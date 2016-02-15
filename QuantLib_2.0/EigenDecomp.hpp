//
//  EigenDecomp.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/13/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef EigenDecomp_hpp
#define EigenDecomp_hpp

#include <stdio.h>
#include <Eigen/Dense>

class EigenDecomp{
// EigenDecomp wraps the Eigen::EigenSolver. It takes input of a real square matrix A, and conducts eigenvalue decomp.
// For a square real matrix A, there is a unique decomp A=VDV^-1. D is diag matrix containing the eigenvalue of A. The column of V is the eigenvector.
// Note: The condition of eigen decomp is: For n*n real matrix A, it must have n linear independent eigen vectors.
    
protected:
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    
public:
    EigenDecomp(){};
    EigenDecomp(const Eigen::MatrixXd &A, bool computeEigenVec= true ) : es(Eigen::EigenSolver<Eigen::MatrixXd>(A, computeEigenVec)){};
    virtual ~EigenDecomp(){};
    
    bool isEigenValueReal() {
        Eigen::MatrixXd EigenValue_img(es.eigenvalues().imag());
        if(EigenValue_img.lpNorm<Eigen::Infinity>()< 1e-8) return true;
        else return false;
    };
    
    Eigen::VectorXd eigenvalues_real(){
        return this->eigenvalues_complex().real();
    };
    
    Eigen::VectorXcd eigenvalues_complex(){
        return es.eigenvalues();
    };
    
    Eigen::MatrixXcd eigenvectors_complex(){
        // return the matrix V
        return es.eigenvectors();
    };
    
    Eigen::MatrixXd eigenvectors_real(){
        // return the real part of V;
        return this->eigenvectors_complex().real();
    };
    

};


#endif /* EigenDecomp_hpp */
