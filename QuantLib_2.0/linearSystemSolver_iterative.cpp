//
//  linearSystemSolver_iterative.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "linearSystemSolver_iterative.hpp"

linearSolver_Jacobi::linearSolver_Jacobi(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations, const Eigen::MatrixXd & GUESS): linearSolver(mtxA, mtxB, tolerance, MaxIterations, GUESS){
    
    
    long nrow= A.cols();
    long ncol= B.cols();
    X.setZero(nrow , ncol);
    if(guess.rows()!=nrow || guess.cols()!= ncol) {guess= X; };
    
    Eigen::MatrixXd E=A.diagonal().asDiagonal().inverse();
    Eigen::MatrixXd R(A); R.diagonal().setZero();
    
    iter_count=0;
    success=true;
    for (long j=0 ; j<ncol; j++){
        long count=0;
        Eigen::VectorXd b(B.col(j));
        Eigen::VectorXd x_old(guess.col(j));
        Eigen::VectorXd x_new;
        double error=1000;
        
        while (error> tol && count< iter_LIMIT) {
            x_new= E* (b-R*x_old);
            error= (x_new- x_old).norm();
            x_old= x_new;
            count++;
            
        }
        X.col(j)= x_new;
        if (count>=iter_LIMIT) {
            success=false;
        }
        iter_count= fmax(iter_count, count);
        
    }
}

linearSolver_WJacobi::linearSolver_WJacobi(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations, double weight, const Eigen::MatrixXd & GUESS): linearSolver(mtxA, mtxB, tolerance, MaxIterations, GUESS){
    
    w= weight;
    long nrow= A.cols();
    long ncol= B.cols();
    X.setZero(nrow , ncol);
    if(guess.rows()!=nrow || guess.cols()!= ncol) guess= X;
    
    Eigen::MatrixXd E(A.diagonal().asDiagonal().inverse());
    Eigen::MatrixXd R(A); R.diagonal().setZero();
    
    iter_count=0;
    success=true;
    for (long j=0 ; j<ncol; j++){
        long count=0;
        Eigen::VectorXd b(B.col(j));
        Eigen::VectorXd x_old(guess.col(j));
        Eigen::VectorXd x_new;
        double error=1000;
        
        while (error> tol && count< iter_LIMIT) {
            x_new= E* (b-R*x_old)*w+ (1-w)*x_old;
            error= (x_new- x_old).norm();
            x_old= x_new;
            count++;
        }
        X.col(j)= x_new;
        if (count>=iter_LIMIT) {
            success=false;
        }
        iter_count= fmax(iter_count, count);
        
    }
    
}

linearSolver_GS::linearSolver_GS(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations, const Eigen::MatrixXd & GUESS): linearSolver(mtxA, mtxB, tolerance, MaxIterations, GUESS){
  
    long nrow= A.cols();
    long ncol= B.cols();
    X.setZero(nrow , ncol);
    if(guess.rows()!=nrow || guess.cols()!= ncol) guess= X;
    
    Eigen::MatrixXd L(A.triangularView<Eigen::Lower>());
    Eigen::MatrixXd U(A.triangularView<Eigen::StrictlyUpper>());
    L= L.inverse();
    
    iter_count=0;
    success=true;
    for (long j=0 ; j<ncol; j++){
        long count=0;
        Eigen::VectorXd b(B.col(j));
        Eigen::VectorXd x_old(guess.col(j));
        Eigen::VectorXd x_new;
        double error=1000;
        
        while (error> tol && count< iter_LIMIT) {
            x_new= L* (b-U*x_old);
            error= (x_new- x_old).norm();
            x_old= x_new;
            count++;
        }
        X.col(j)= x_new;
        if (count>=iter_LIMIT) {
            success=false;
        }
        iter_count= fmax(iter_count, count);
        
    }
    
    
}


linearSolver_SOR::linearSolver_SOR(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd &mtxB, double tolerance, long MaxIterations, double weight, const Eigen::MatrixXd & GUESS): linearSolver(mtxA, mtxB, tolerance, MaxIterations, GUESS) {
  
    w= weight;
    long nrow= A.cols();
    long ncol= B.cols();
    X.setZero(nrow , ncol);
    if(guess.rows()!=nrow || guess.cols()!= ncol) guess= X;
    
    Eigen::MatrixXd D(A.diagonal().asDiagonal());
    Eigen::MatrixXd L(A.triangularView<Eigen::StrictlyLower>());
    Eigen::MatrixXd U(A.triangularView<Eigen::StrictlyUpper>());
    
    Eigen::MatrixXd S=(D+L*w).inverse();
    Eigen::MatrixXd R= w*U+ (w-1)*D;
    
    
    iter_count=0;
    success=true;
    for (long j=0 ; j<ncol; j++){
        long count=0;
        Eigen::VectorXd b(B.col(j));
        Eigen::VectorXd x_old(guess.col(j));
        Eigen::VectorXd x_new;
        double error=1000;
        
        while (error> tol && count< iter_LIMIT) {
            x_new= S* (w*b-R*x_old);
            error= (x_new- x_old).norm();
            x_old= x_new;
            count++;
        }
        X.col(j)= x_new;
        if (count>=iter_LIMIT) {
            success=false;
        }
        iter_count= fmax(iter_count, count);
        
    }

}
