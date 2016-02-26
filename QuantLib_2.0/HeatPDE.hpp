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
    
    
    virtual void ForwardEuler();
    virtual void BackwardEuler();
    virtual void CrankNicolson();
    // backwardEuler and CrankNicolson involves solving linear system. We also implement a version using iterative linear system solver (SOR)
    virtual void BackwardEuler_SOR(double weight, double tol, long iterMAX);
    virtual void CrankNicolson_SOR(double weight, double tol, long iterMAX);
    
    
public:
    //constructors
    HeatPDE(){};
    HeatPDE(double coef, double xmin, double xmax, double tmax,
            long x_steps, long t_steps,
            const ROOT::Math::IBaseFunctionOneDim & down,
            const ROOT::Math::IBaseFunctionOneDim & left,
            const ROOT::Math::IBaseFunctionOneDim & right,
            int algorithm, double weight= 1.1, double tol= 1.e-5, long iterMAX= 50);

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
    virtual void ForwardEuler();
    virtual void BackwardEuler();
    virtual void CrankNicolson();
    // backwardEuler and CrankNicolson involves solving linear system. We also implement a version using iterative linear system solver (SOR)
    virtual void BackwardEuler_SOR(double weight, double tol, long iterMAX);
    virtual void CrankNicolson_SOR(double weight, double tol, long iterMAX);
public:
    // constructors
    HeatPDE_gridUpdate(){};
    HeatPDE_gridUpdate(double coef, double xmin, double xmax, double tmax,
                       long x_steps, long t_steps,
                       const ROOT::Math::IBaseFunctionOneDim & down,
                       const ROOT::Math::IBaseFunctionOneDim & left,
                       const ROOT::Math::IBaseFunctionOneDim & right,
                       const ROOT::Math::IBaseFunctionMultiDim & lowbound,
                       int algorithm, double weight= 1.1, double tol= 1.e-5, long iterMAX= 50);
    // destructors
    virtual ~ HeatPDE_gridUpdate(){};
    
};

#endif /* HeatPDE_hpp */
