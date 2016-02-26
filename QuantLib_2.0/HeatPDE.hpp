//
//  HeatPDE.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/20/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef HeatPDE_hpp
#define HeatPDE_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <Math/IFunction.h>
#include <math.h>
#include <Math/Integrator.h>
#include "LU_FullPiv.hpp"
#include "linearSystemSolver_iterative.hpp"


namespace QLib {
    namespace FD{
        
        using namespace QLib::FD;
        using namespace QLib;
class HeatPDE{
// HeatPDE defines the interface and methods of the numerical solutions of one dimesional standard heat PDE problem.
protected:
    //Standard heat PDE is in the form of
    // u_t= k*u_xx, with boundary condition u(x,0)= f(x), x in R and t>0. We are interested in u(x,T) for some T>0.
    // For the numerical solution we need to constrain the expansion of x besides the boundary of t=0, so another boundary function g_left and g_right is defined.
    double k;
    const ROOT::Math::IBaseFunctionOneDim * f;
    const ROOT::Math::IBaseFunctionOneDim * g_left;
    const ROOT::Math::IBaseFunctionOneDim * g_right;

    
    // The numerical solution conducts grid computation. We sometimes interested in the grid value, and use a matrix to store the computation.
    long N, M; // N decides how many intervals we discretize the x domain, and M decides the t domain.
    double xMIN, xMAX, tMAX; // xMIN, xMAX and T defines the solution domain.
    Eigen::MatrixXd grid; // grid is a matrix of (M+1)*(N+1). The first row corresponds to the time 0. The first column corresponds to the xMIN. 
    
    // There are three popular methods for the grid computation, Forward Euler, Backward Euler and Crank-Nicosan. We code it as 0, 1 and 2. Code 3 means BackwardEuler_SOR. Code 4 means CrankNicolson_SOR
    int method;
    double alpha; // alpha is a critical param deciding how good the convergence is.
    bool status; // Backward and CrankNicolson involves linear system solution. status tracks whether there is any unexpected computation error. status= true if everything good, flase if else.
    
    
    virtual void ForwardEuler(){
        for (long i=1; i<M+1; i++) {
            for (long j=1; j<N; j++) {
                grid(i,j)= alpha* (grid(i-1, j-1)+ grid(i-1, j+1)) + (1-2*alpha)* grid(i-1, j);
            }
        }
    };
    virtual void BackwardEuler(){
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+2*alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        LU_FullPiv lu_f(A);
        Eigen::PermutationMatrix<Eigen::Dynamic> P= lu_f.matrixP();
        Eigen::PermutationMatrix<Eigen::Dynamic> Q= lu_f.matrixQ();
        Eigen::MatrixXd L= lu_f.matrixL();
        Eigen::MatrixXd U= lu_f.matrixU();
        
        for(long i=1; i<M+1; i++){
            Eigen::VectorXd b= grid.block(i-1, 1, 1, N-1).transpose();
            b(0)+= alpha* grid(i, 0);
            b(N-2)+= alpha* grid(i, N);
            
            Eigen::VectorXd y= L.triangularView<Eigen::Lower>().solve(P*b);
            Eigen::VectorXd x= U.triangularView<Eigen::Upper>().solve(y);
            if ((y-U.triangularView<Eigen::Upper>()* x).lpNorm<Eigen::Infinity>()> 1e-6) status= false;
            
            x= Q*x;
            
            grid.block(i, 1, 1, N-1)= x.transpose();
        }

    };
    virtual void CrankNicolson() {
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha/2.0;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        Eigen::MatrixXd B; B.setZero(N-1, N-1);
        diag.setOnes(N-1);
        diag*= (1-alpha) ;
        B.diagonal()= diag;
        diag.setOnes(N-2); diag*=alpha/2.0;
        B.diagonal(1)=diag;
        B.diagonal(-1)=diag;
        
        LU_FullPiv lu_f(A);
        Eigen::PermutationMatrix<Eigen::Dynamic> P= lu_f.matrixP();
        Eigen::PermutationMatrix<Eigen::Dynamic> Q= lu_f.matrixQ();
        Eigen::MatrixXd L= lu_f.matrixL();
        Eigen::MatrixXd U= lu_f.matrixU();
        
        for (long i=1; i<M+1; i++) {
            Eigen::VectorXd b=grid.block(i-1, 1, 1, N-1).transpose();
            b= B* b;
            
            b(0)+=      alpha/2 * (grid(i,0)+ grid(i-1,0));
            b(N-2)+=    alpha/2 * (grid(i,N)+ grid(i-1,N));
            
            Eigen::VectorXd y= L.triangularView<Eigen::Lower>().solve(P*b);
            Eigen::VectorXd x= U.triangularView<Eigen::Upper>().solve(y);
            if ((y-U.triangularView<Eigen::Upper>()* x).lpNorm<Eigen::Infinity>()> 1.0e-6) status= false;
            
            x= Q*x;
            
            grid.block(i, 1, 1, N-1)= x.transpose();
            
        }
        

    };
    // backwardEuler and CrankNicolson involves solving linear system. We also implement a version using iterative linear system solver (SOR)
    
    virtual void BackwardEuler_SOR(double weight, double tol, long iterMAX){
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+2*alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        for(long i=1; i<M+1; i++){
            Eigen::VectorXd b= grid.block(i-1, 1, 1, N-1).transpose();
            b(0)+= alpha* grid(i, 0);
            b(N-2)+= alpha* grid(i, N);
            
            linearSolver_SOR l_sor(A, b, tol, iterMAX, weight);
            grid.block(i, 1, 1, N-1)=  l_sor.solve().transpose();
            status= l_sor.convergence();
            
            std::cout<< i<< std::endl;
        }
        
        
    };
    
    virtual void CrankNicolson_SOR(double weight, double tol, long iterMAX) {
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha/2.0;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        Eigen::MatrixXd B; B.setZero(N-1, N-1);
        diag.setOnes(N-1);
        diag*= (1-alpha) ;
        B.diagonal()= diag;
        diag.setOnes(N-2); diag*=alpha/2.0;
        B.diagonal(1)=diag;
        B.diagonal(-1)=diag;
        
        for (long i=1; i<M+1; i++) {
            Eigen::VectorXd b=grid.block(i-1, 1, 1, N-1).transpose();
            b= B* b;
            
            b(0)+=      alpha/2 * (grid(i,0)+ grid(i-1,0));
            b(N-2)+=    alpha/2 * (grid(i,N)+ grid(i-1,N));
            
            linearSolver_SOR l_sor(A, B, tol, iterMAX, weight);
            
            grid.block(i, 1, 1, N-1)= l_sor.solve().transpose();
            
            status= l_sor.convergence();
        }

    };
    
    
public:
    //constructors
    HeatPDE(){};
    HeatPDE(double coef, double xmin, double xmax, double tmax,
            long x_steps, long t_steps,
            const ROOT::Math::IBaseFunctionOneDim & down,
            const ROOT::Math::IBaseFunctionOneDim & left,
            const ROOT::Math::IBaseFunctionOneDim & right,
            int algorithm, double weight= 1.1, double tol= 1.e-5, long iterMAX= 50){
        k= coef;
        xMIN= xmin;
        xMAX= xmax;
        tMAX= tmax;
        N= x_steps;
        M= t_steps;
        f= &down;
        g_left= &left;
        g_right=&right;
        method= algorithm;
        
        double dt= tMAX/ M;
        double dx= fabs(xMAX- xMIN)/ N;
        
        alpha= dt/(dx* dx);
        grid.setZero(M+1, N+1);
        // compute first row of grid
        for (long j=0; j<N+1; j++) {
            grid(0,j)= f->operator()(xMIN+ dx* j);
        }
        // compute first and last column of grid
        for (long i=0; i<M+1; i++) {
            grid(i,0)= g_left-> operator()(dt*i);
            grid(i,N)= g_right-> operator()(dt*i);
        }
        
        status= true;
        
        if (method==0){
            this->ForwardEuler();
        }else if( method==1){
            this->BackwardEuler();
        }else if (method==2){
            this->CrankNicolson();
        }else if(method==3){
            this->BackwardEuler_SOR(weight, tol, iterMAX);
        }else if(method==4){
            this->CrankNicolson_SOR(weight, tol, iterMAX);
        }else {
            std::cout<< "Heat PDE solver: undefined algorithm.\n\n";
        }
        
    };

    //destructor
    virtual ~ HeatPDE(){};
    
    //methods
    virtual double dx(){return fabs(xMAX- xMIN)/N;};
    virtual double dt(){return tMAX/M; };
    virtual double heatPDECoef() const {return k;};
    virtual Eigen::MatrixXd computationGrid() const {return grid;};
    virtual int computationMethod() const  {return method;};
    virtual double computationAlpha() const {return alpha;};
    virtual bool computationSuccess() const {return status;};
    virtual double evaluation(double x, double t) const {
        // evaluate the result at (x,t). x,t must be in [xMIN, xMAX]*[0,tMAX].
        
        double dx= fabs(xMAX- xMIN)/ N;
        double dt= tMAX/ M;
        
        long j= (x-xMIN)/ dx;
        long i= t/dt;
        
        double w11, w12, w21, w22;
        w11= (t- i*dt)/ dt;
        w12= w11* (x- xMIN- dx* j)/dx;
        w11= w11- w12;
        w21= 1- w11 - w12;
        w22= w21* (x- xMIN- dx* j)/dx;
        w21= w21- w22;
        
        if (j+1> N){
            if(i+1>M) {
                return grid(i,j);
            }else
                return (w11+w12)* grid(i+1,j)+ (w21+w22)* grid(i,j);
        }else if (i+1 >M) {
            if (j+1>N){
                return grid(i,j);
            }else
                return (w11+w21)*grid(i,j)+ (w12+w22)*grid(i,j+1);
        }else {
            double tmp= grid(i,j)* w21+ w22* grid(i,j+1)+ w11* grid(i+1,j)+ w12* grid(i+1, j+1);
            return tmp;
        }
    };
    
};






class HeatPDE_gridUpdate: public HeatPDE {
    // HeatPDE_gridUpdate accommodates the American option pricing PDE kernals. It introduce a lower bounds function which dynamically sets the lower bounds of u(x,t) at each computation nodes
    const ROOT::Math::IBaseFunctionMultiDim * lowBound;
    
    // methods
    virtual void ForwardEuler(){
        double dx= fabs(xMAX- xMIN)/ N;
        double dt= tMAX/ M;
        
        for (long i=1; i<M+1; i++) {
            for (long j=1; j<N; j++) {
                grid(i,j)= alpha* (grid(i-1, j-1)+ grid(i-1, j+1)) + (1-2*alpha)* grid(i-1, j);
                double y[2]={xMIN+ j*dx,  dt*i};
                grid(i,j)= fmax(grid(i,j), lowBound->operator()(y));
            }
        }

    };
    virtual void BackwardEuler(){
        double dx= fabs(xMAX- xMIN)/ N;
        double dt= tMAX/ M;
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+2*alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        LU_FullPiv lu_f(A);
        Eigen::PermutationMatrix<Eigen::Dynamic> P= lu_f.matrixP();
        Eigen::PermutationMatrix<Eigen::Dynamic> Q= lu_f.matrixQ();
        Eigen::MatrixXd L= lu_f.matrixL();
        Eigen::MatrixXd U= lu_f.matrixU();
        
        for(long i=1; i<M+1; i++){
            Eigen::VectorXd b= grid.block(i-1, 1, 1, N-1).transpose();
            b(0)+= alpha* grid(i, 0);
            b(N-2)+= alpha* grid(i, N);
            
            Eigen::VectorXd y= L.triangularView<Eigen::Lower>().solve(P*b);
            Eigen::VectorXd x= U.triangularView<Eigen::Upper>().solve(y);
            if ((y-U.triangularView<Eigen::Upper>()* x).lpNorm<Eigen::Infinity>()> 1e-6) status= false;
            
            x= Q*x;
            
            grid.block(i, 1, 1, N-1)= x.transpose();
            for (long j=1; j<N; j++) {
                double y[2]= {xMIN+ dx* j, dt*i};
                grid(i,j)= fmax(grid(i,j), lowBound-> operator()(y));
            }
        }
        

    };
    virtual void CrankNicolson(){
        double dx= fabs(xMAX- xMIN)/ N;
        double dt= tMAX/ M;
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha/2.0;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        Eigen::MatrixXd B; B.setZero(N-1, N-1);
        diag.setOnes(N-1);
        diag*= (1-alpha) ;
        B.diagonal()= diag;
        diag.setOnes(N-2); diag*=alpha/2.0;
        B.diagonal(1)=diag;
        B.diagonal(-1)=diag;
        
        LU_FullPiv lu_f(A);
        Eigen::PermutationMatrix<Eigen::Dynamic> P= lu_f.matrixP();
        Eigen::PermutationMatrix<Eigen::Dynamic> Q= lu_f.matrixQ();
        Eigen::MatrixXd L= lu_f.matrixL();
        Eigen::MatrixXd U= lu_f.matrixU();
        
        for (long i=1; i<M+1; i++) {
            Eigen::VectorXd b=grid.block(i-1, 1, 1, N-1).transpose();
            b= B* b;
            
            b(0)+=      alpha/2 * (grid(i,0)+ grid(i-1,0));
            b(N-2)+=    alpha/2 * (grid(i,N)+ grid(i-1,N));
            
            Eigen::VectorXd y= L.triangularView<Eigen::Lower>().solve(P*b);
            Eigen::VectorXd x= U.triangularView<Eigen::Upper>().solve(y);
            if ((y-U.triangularView<Eigen::Upper>()* x).lpNorm<Eigen::Infinity>()> 1.0e-6) status= false;
            
            x= Q*x;
            
            grid.block(i, 1, 1, N-1)= x.transpose();
            
            for (long j=1; j<N; j++) {
                double y[2]= {xMIN+ dx* j, dt*i};
                grid(i,j)= fmax(grid(i,j), lowBound-> operator()(y));
            }
        }

    };
    // backwardEuler and CrankNicolson involves solving linear system. We also implement a version using iterative linear system solver (SOR)
    virtual void BackwardEuler_SOR(double weight, double tol, long iterMAX) {
        double dx= fabs(xMAX- xMIN)/N;
        double dt= tMAX/M;
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+2*alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        for(long i=1; i<M+1; i++){
            Eigen::VectorXd b= grid.block(i-1, 1, 1, N-1).transpose();
            b(0)+= alpha* grid(i, 0);
            b(N-2)+= alpha* grid(i, N);
            
            linearSolver_SOR l_sor(A, b, tol, iterMAX, weight);
            grid.block(i, 1, 1, N-1)=  l_sor.solve().transpose();
            status= l_sor.convergence();
            
            for (long j=1; j<N; j++) {
                double y[2]= {xMIN+ dx* j, dt*i};
                grid(i,j)= fmax(grid(i,j), lowBound-> operator()(y));
            }
        }

    };
    virtual void CrankNicolson_SOR(double weight, double tol, long iterMAX){
        double dx= fabs(xMAX- xMIN)/N;
        double dt= tMAX/M;
        
        Eigen::MatrixXd A; A.setZero(N-1, N-1);
        Eigen::VectorXd diag; diag.setOnes(N-1);
        diag*= (1+alpha) ;
        A.diagonal()= diag;
        diag.setOnes(N-2); diag*=-alpha/2.0;
        A.diagonal(1)=diag;
        A.diagonal(-1)=diag;
        
        Eigen::MatrixXd B; B.setZero(N-1, N-1);
        diag.setOnes(N-1);
        diag*= (1-alpha) ;
        B.diagonal()= diag;
        diag.setOnes(N-2); diag*=alpha/2.0;
        B.diagonal(1)=diag;
        B.diagonal(-1)=diag;
        
        for (long i=1; i<M+1; i++) {
            Eigen::VectorXd b=grid.block(i-1, 1, 1, N-1).transpose();
            b= B* b;
            
            b(0)+=      alpha/2 * (grid(i,0)+ grid(i-1,0));
            b(N-2)+=    alpha/2 * (grid(i,N)+ grid(i-1,N));
            
            linearSolver_SOR l_sor(A, B, tol, iterMAX, weight);
            
            grid.block(i, 1, 1, N-1)= l_sor.solve().transpose();
            
            status= l_sor.convergence();
            for (long j=1; j<N; j++) {
                double y[2]= {xMIN+ dx* j, dt*i};
                grid(i,j)= fmax(grid(i,j), lowBound-> operator()(y));
            }
        }

    };
public:
    // constructors
    HeatPDE_gridUpdate(){};
    HeatPDE_gridUpdate(double coef, double xmin, double xmax, double tmax,
                       long x_steps, long t_steps,
                       const ROOT::Math::IBaseFunctionOneDim & down,
                       const ROOT::Math::IBaseFunctionOneDim & left,
                       const ROOT::Math::IBaseFunctionOneDim & right,
                       const ROOT::Math::IBaseFunctionMultiDim & lowbound,
                       int algorithm, double weight= 1.1, double tol= 1.e-5, long iterMAX= 50){
        k= coef;
        xMIN= xmin;
        xMAX= xmax;
        tMAX= tmax;
        N= x_steps;
        M= t_steps;
        f= &down;
        g_left= &left;
        g_right=&right;
        lowBound= & lowbound;
        method= algorithm;
        
        double dt= tMAX/ M;
        double dx= fabs(xMAX- xMIN)/ N;
        
        alpha= dt/(dx* dx);
        grid.setZero(M+1, N+1);
        // compute first row of grid
        for (long j=0; j<N+1; j++) {
            grid(0,j)= f->operator()(xMIN+ dx* j);
        }
        // compute first and last column of grid
        for (long i=0; i<M+1; i++) {
            grid(i,0)= g_left-> operator()(dt*i);
            grid(i,N)= g_right-> operator()(dt*i);
        }
        
        status= true;
        
        if (method==0){
            this->ForwardEuler();
        }else if( method==1){
            this->BackwardEuler();
        }else if (method==2){
            this->CrankNicolson();
        }else if(method==3){
            this->BackwardEuler_SOR(weight, tol, iterMAX);
        }else if(method==4){
            this->CrankNicolson_SOR(weight, tol, iterMAX);
        }else {
            std::cout<< "Heat PDE solver: undefined algorithm.\n\n";
        }
        

    };
    // destructors
    virtual ~ HeatPDE_gridUpdate(){};
    
};
    }
}

#endif /* HeatPDE_hpp */
