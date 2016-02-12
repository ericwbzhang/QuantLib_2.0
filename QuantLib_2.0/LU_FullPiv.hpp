//
//  LU_FullPiv.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef LU_FullPiv_hpp
#define LU_FullPiv_hpp

#include <stdio.h>
#include <Eigen/Dense>

class LU_FullPiv{
// LU_FullPiv wraps the Eigen::FullPivLU class. It gets initialized with an input matrix A to be decomposed. A is a general real matix, not necessary to be square matrix.
// For a general matrix A with size m*n, it can be decomposed as PAQ= LU, P and Q are permutation matrix. L is a m*m unit-lower triangular matrix, while U is an upper triangular matrix with size m*n.

protected:
    Eigen::FullPivLU<Eigen::MatrixXd > lu_f;
    
public:
    LU_FullPiv(){};
    LU_FullPiv(const Eigen::MatrixXd & A ) {
        lu_f= Eigen::FullPivLU<Eigen::MatrixXd> (A);
    }
    
    virtual ~LU_FullPiv(){};
  
public:
    
    Eigen::PermutationMatrix<Eigen::Dynamic> matrixP() {return lu_f.permutationP();};
    Eigen::PermutationMatrix<Eigen::Dynamic> matrixQ() {return lu_f.permutationQ();};
    Eigen::MatrixXd matrixL(){
        Eigen::MatrixXd LU= lu_f.matrixLU();
        long nrow= LU.rows();
        long ncol= LU.cols();
        Eigen::MatrixXd L; L.setIdentity(nrow, nrow);
    
        L.block(0,0, nrow, nrow<ncol? nrow: ncol).triangularView<Eigen::StrictlyLower>()=  LU.block(0,0,nrow, fmin(nrow, ncol));
        return L;
    };
    Eigen::MatrixXd matrixU(){
        Eigen::MatrixXd LU= lu_f.matrixLU();
        long nrow= LU.rows();
        long ncol= LU.cols();
        Eigen::MatrixXd U(LU.triangularView<Eigen::Upper>());
        U.conservativeResize(nrow, ncol); 
        return U;
    };
    
    bool isInvertible(){return lu_f.isInvertible();};
    double determinant(){return lu_f.determinant();};
    
};


#endif /* LU_FullPiv_hpp */
