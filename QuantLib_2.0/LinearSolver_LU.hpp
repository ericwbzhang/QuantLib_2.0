//
//  LinearSolver_LU.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/21/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef LinearSolver_LU_hpp
#define LinearSolver_LU_hpp

#include <stdio.h>
#include <Eigen/Dense>

class LinearSolver_LU{
// LinearSolver_LU defines the interface and method to solve a linear system L*U* x= B, where L is a lower m*m triangular matrix and U is an upper m*n triangular matrix.
// Note: If L&U comes from a LU decomposition, L is invertible and has size m*m with diagonal 1, U is a m*n matrix.
    
protected:
    
public:
    LinearSolver_LU(){};
//    LinearSolver_LU(const Eigen::MatrixXd & L, const Eigen::MatrixXd & U){
//    
//    };
//    
    virtual ~ LinearSolver_LU(){};
    
    bool solve(const Eigen::MatrixXd & L, const Eigen::MatrixXd &U, const Eigen::MatrixXd &B,  Eigen::MatrixXd &X){
        // solve takes the linear system as input and update the matrix X. It returns true if the system has a solution and false if the system has no solution (X min the L2 norm of B-LUX)
        Eigen::MatrixXd Y= L.triangularView<Eigen::Lower>().solve(B);
        X= U.triangularView<Eigen::Upper>().solve(Y);
        if ((Y- U.triangularView<Eigen::Upper>()* X).lpNorm<Eigen::Infinity>()< 1e-8) return true;
        else return false;
    };
};

#endif /* LinearSolver_LU_hpp */
