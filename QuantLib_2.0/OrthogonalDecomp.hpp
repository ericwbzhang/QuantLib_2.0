//
//  OrthogonalDecomp.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/13/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef OrthogonalDecomp_hpp
#define OrthogonalDecomp_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include "EigenDecomp.hpp"

class OrthogonalDecomp{
// OrthogonalDecomp conducts Orthogonal decomp for a real symmetric matrix A. A= VDV'. V is orthogonal matrix. D is diagonal.
// It takes use of Eigen::HouseholderQR method

protected:
    bool symmetric;
    Eigen::MatrixXd V;
    Eigen::MatrixXd D;
    
public:
    OrthogonalDecomp(){};
    OrthogonalDecomp(const Eigen::MatrixXd &A){
        symmetric= true;
        if ((A-A.transpose()).lpNorm<Eigen::Infinity>() >1e-8) {
            symmetric=false;
        }
        
        if(symmetric){
            EigenDecomp ed(A);
            D= ed.eigenvalues_real().asDiagonal();
            V= ed.eigenvectors_real();
            V= V.householderQr().householderQ();
        }
    };
    virtual ~OrthogonalDecomp(){};
    
    bool isSymmetric(){return symmetric;};
    
    Eigen::MatrixXd matrixV(){return V;};
    Eigen::MatrixXd MatrixD(){return D;};

};
#endif /* OrthogonalDecomp_hpp */
