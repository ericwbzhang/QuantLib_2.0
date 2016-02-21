//
//  HeatPDE.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/20/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "HeatPDE.hpp"
#include "LU_FullPiv.hpp"
#include "linearSystemSolver_iterative.hpp"

HeatPDE:: HeatPDE(double coef, double xmin, double xmax, double tmax,
                  long x_steps, long t_steps,
                  const ROOT::Math::IBaseFunctionOneDim & down,
                  const ROOT::Math::IBaseFunctionOneDim & left,
                  const ROOT::Math::IBaseFunctionOneDim & right,
                  int algorithm, double weight, double tol, long iterMAX){
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

}


void HeatPDE::ForwardEuler(){
    double tmp;
    for (long i=1; i<M+1; i++) {
        for (long j=1; j<N; j++) {
            grid(i,j)= alpha* (grid(i-1, j-1)+ grid(i-1, j+1)) + (1-2*alpha)* grid(i-1, j);
            tmp= grid(i,j);
        }
    }
}


void HeatPDE::BackwardEuler(){
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
}



void HeatPDE::CrankNicolson(){
    
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

}

void HeatPDE::BackwardEuler_SOR(double weight, double tol, long iterMAX){
    
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
    
    
}



void HeatPDE::CrankNicolson_SOR(double weight, double tol, long iterMAX){
    
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
}