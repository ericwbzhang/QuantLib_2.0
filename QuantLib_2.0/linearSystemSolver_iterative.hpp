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
    linearSolver_Jacobi(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations=50, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd());

};

// Weighted Jacobi Solver

class linearSolver_WJacobi: public linearSolver{
    // Note: for Weighted-Jacobi method, convergence depends on matrix A and w. If A is irreducible and diagonally dominant, Jacobi will converge, and we usually set w=2/3.0.
    // iter_count counts maximum of iterations (usually no greater than 20) Jacobi method spends for each column of B. iter_LIMIT set the upper limit of iter_count.
protected:
    double w;
    
public:
    linearSolver_WJacobi(){};
    linearSolver_WJacobi(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd &mtxB, double tolerance, long MaxIterations=50, double weight= 2.0/3.0, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd());
    virtual ~linearSolver_WJacobi(){};
};


// Gauss-Seidel

class linearSolver_GS: public linearSolver{
    // Mote: for GS method, convergence depends on A. If A is either spd or irreducibly diagonally dominant, GS converges.
    // iter_count tracks the max of iteratiions GS method spends for each column of B. iter_LIMIT set the upper limit.
    
public:
    linearSolver_GS(){};
    virtual ~linearSolver_GS(){};
    linearSolver_GS(const Eigen::MatrixXd & mtxA, const Eigen::MatrixXd & mtxB, double tolerance, long MaxIterations=50, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd());
    
};


// SOR

class linearSolver_SOR: public linearSolver{
// Note: for SOR the convergence depends on A and w. If A is spd and w is in (0,2), SOR converges. Note when w=1, SOR is the same with GS.
    
protected:
    double w;
public:
    linearSolver_SOR(){};
    virtual ~linearSolver_SOR(){};
    linearSolver_SOR(const Eigen::MatrixXd &mtxA, const Eigen::MatrixXd &mtxB, double tolerance, long MaxIterations=50, double weight= 1.1, const Eigen::MatrixXd & GUESS= Eigen::MatrixXd());

};

#endif /* linearSystemSolver_iterative_hpp */
