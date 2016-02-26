//
//  LU_ParPiv.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef LU_ParPiv_hpp
#define LU_ParPiv_hpp

#include <stdio.h>
#include <Eigen/Dense>
namespace QLib{
    using namespace QLib;

class LU_ParPiv{
// LU_ParPiv wraps Eigen::PartialPivLU. It plays partial permutation LU decomp to the input matrix A, which must be an invertible square matrix.
// For an invertible square matrix A, there is unique decomp PA= LU, where P is a permutation matrix and L is unit-lower triangular matrix and U is an upper triangular matrix.
    
protected:
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_p;
    
public:
    LU_ParPiv(){};
    LU_ParPiv(const Eigen::MatrixXd &A): lu_p(A.lu()){};
    
    virtual ~LU_ParPiv(){};
    
    Eigen::PermutationMatrix<Eigen::Dynamic> matrixP(){return lu_p.permutationP();};
    Eigen::MatrixXd matrixL () {
        Eigen::MatrixXd LU= lu_p.matrixLU();
        Eigen::MatrixXd L; L.setIdentity(LU.rows(), LU.cols());
        L.triangularView<Eigen::StrictlyLower>()= LU;
        return L;
    };
    
    Eigen::MatrixXd matrixU(){
        Eigen::MatrixXd LU= lu_p.matrixLU();
        Eigen::MatrixXd U( LU.triangularView<Eigen::Upper>());
        return U;
    };
};

}


#endif /* LU_ParPiv_hpp */
