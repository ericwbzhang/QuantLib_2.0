//
//  linearSystemSolver_iterative.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef linearSystemSolver_iterative_hpp
#define linearSystemSolver_iterative_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <iostream>


namespace QLib {
    using namespace QLib;

// Base Class
class linearSolver {

protected:
    long double tol;
    long iter_count;
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd X;
    
    Eigen::MatrixXd guess;
    long iter_LIMIT;
    
    bool success;
    // succeed to convergence or not
public:
    
    linearSolver(){};
    linearSolver(const Eigen::MatrixXd &mtxA, const  Eigen::MatrixXd &mtxB,  double tolerance, long MaxIteratiions=50, const Eigen::MatrixXd & GUESS=Eigen::MatrixXd()):A(mtxA), B(mtxB), tol(tolerance), iter_LIMIT(MaxIteratiions), guess(GUESS){};
    virtual ~linearSolver(){};
    
    bool convergence(){return success;};
    long iterations(){return iter_count;};
    Eigen::MatrixXd solve(){return X;};
    long iter_limit() {return iter_LIMIT;};

};


// Jacobi Solver
class linearSolver_Jacobi : public linearSolver{
// Note: for Jacobi method, convergence depends on matrix A. If A is irreducible and diagonally dominant, Jacobi will converge.
// iter_count counts maximum of iterations (usually no greater than 20) Jacobi method spends for each column of B. iter_LIMIT set the upper limit of iter_count.

public:
    virtual ~linearSolver_Jacobi(){};
    linearSolver_Jacobi(){};
    linearSolver_Jacobi(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations=50, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd()){
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

    };

};

// Weighted Jacobi Solver

class linearSolver_WJacobi: public linearSolver{
    // Note: for Weighted-Jacobi method, convergence depends on matrix A and w. If A is irreducible and diagonally dominant, Jacobi will converge, and we usually set w=2/3.0.
    // iter_count counts maximum of iterations (usually no greater than 20) Jacobi method spends for each column of B. iter_LIMIT set the upper limit of iter_count.
protected:
    double w;
    
public:
    linearSolver_WJacobi(){};
    linearSolver_WJacobi(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd &mtxB, double tolerance, long MaxIterations=50, double weight= 2.0/3.0, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd()){
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

    };
    virtual ~linearSolver_WJacobi(){};
};


// Gauss-Seidel

class linearSolver_GS: public linearSolver{
    // Mote: for GS method, convergence depends on A. If A is either spd or irreducibly diagonally dominant, GS converges.
    // iter_count tracks the max of iteratiions GS method spends for each column of B. iter_LIMIT set the upper limit.
    
public:
    linearSolver_GS(){};
    virtual ~linearSolver_GS(){};
    linearSolver_GS(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations=50, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd()){
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

    };
    
};


// SOR

class linearSolver_SOR: public linearSolver{
// Note: for SOR the convergence depends on A and w. If A is spd and w is in (0,2), SOR converges. Note when w=1, SOR is the same with GS.
    
protected:
    double w;
public:
    linearSolver_SOR(){};
    virtual ~linearSolver_SOR(){};
    linearSolver_SOR(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd &mtxB, double tolerance, long MaxIterations=50, double weight= 1.1, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd()){
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
        

    };

};
}
#endif /* linearSystemSolver_iterative_hpp */
