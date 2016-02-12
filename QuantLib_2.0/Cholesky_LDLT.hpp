//
//  Cholesky_LDLT.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/11/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef Cholesky_LDLT_hpp
#define Cholesky_LDLT_hpp

#include <stdio.h>
#include <Eigen/Dense>


class Cholesky_LDLT {

    // This class wraps the Eigen::LDLT class, which performs a robust Cholesky decomp of a symmetric posive semidefinite matrix A, such that A= P'U'DUP where P is a permutation matrix and U is a unit upper triangle matrix and D is a diag matrix
    // The LDLT decomp can also be applied on symmetric negative semidefinite matrix, then A= P'U* D U P, U* is the conjugate tranpose of U.
    // Note the permuation matrix is orthognal by definition, P'P= I
    
protected:
    Eigen::LDLT<Eigen::MatrixXd> ldlt;
    bool spd_flag;
    long n;
    // spd_flag= true if D has positive real entries. It means A is symmetric positive definite.
    
public:
    Cholesky_LDLT(){};
    Cholesky_LDLT(const Eigen::MatrixXd &A) {
        n= A.rows();
        ldlt= A.ldlt();
        Eigen::VectorXd diagD(ldlt.vectorD());
        spd_flag= true;
        for (int i=0; i< diagD.size(); i++) {
            if (diagD(i)<= 0) spd_flag= false;
        }
    };
    
    virtual ~Cholesky_LDLT(){};
    
    Eigen::MatrixXd matrixU() {return ldlt. matrixU();};
    Eigen::MatrixXd matrixD() {
        Eigen::MatrixXd tmp;
        tmp.setIdentity(n,n);
        tmp.diagonal()= ldlt.vectorD();
        return tmp;
    };
    
    Eigen::MatrixXd matrixP() {
        return ldlt.transpositionsP()* Eigen::MatrixXd::Identity(n,n);
    };
    
    bool spd() {return spd_flag;};
    
    Eigen::MatrixXd factorU(){
        Eigen::VectorXd diag_D( ldlt.vectorD()) ;
        for (long i=0; i< diag_D.size(); i++){
            diag_D(i)= sqrt(diag_D(i));
        }
        Eigen::MatrixXd tmp;
        tmp.setIdentity(n,n);
        tmp.diagonal()= diag_D;
        return tmp* this-> matrixU()* this-> matrixP();
    };
    
};

#endif /* Cholesky_LDLT_hpp */
