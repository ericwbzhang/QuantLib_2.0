//
//  SVD_Jb.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef SVD_Jb_hpp
#define SVD_Jb_hpp

#include <stdio.h>
#include <Eigen/Dense>

class SVD_Jb{
// SVD_Jb wraps the Eigen::JacobiSVD class. It takes input  a general real matrix A of n*p. The SVD of A= USV'. U n*n orthogonal matrix, V p*p orthogonal matrix and S is n*p diagonal matrix with non negative entries.
// Note: In case of a rectangular n-by-p matrix, letting m be the smaller value among n and p, there are only m singular vectors; the remaining columns of U and V do not correspond to actual singular vectors. Asking for thin U or V means asking for only their m first columns to be formed. So U is then a n-by-m matrix, and V is then a p-by-m matrix. Notice that thin U and V are all you need for (least squares) solving.
// ComputeFull is the slowest, while ComputeThin can significantly improve the speed especially when the size of A are extremely unbalanced. Not compute the U and V are definitely fastest. 

protected:
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    int fullU; // 1 means we are going to compute the full matrix U, 0 means computing thin matrix U, and -1 means matrix U is not computed.
    int fullV; // 1 means we are going to compute the full matrix V, 0 means computing thin matrix V, and -1 means matrix V is not computed.
    long nrow, ncol;
public:
    SVD_Jb(){};
    SVD_Jb(const Eigen::MatrixXd &A, Eigen::DecompositionOptions computeOptionU, Eigen::DecompositionOptions computeOptionV) {
        nrow= A.rows();
        ncol= A.cols();
        svd= Eigen::JacobiSVD<Eigen::MatrixXd> (A, computeOptionU | computeOptionV);
        fullU=false;
        fullV=false;
        
        if (computeOptionU== Eigen::ComputeFullU) {
            fullU=true;
        }
        if (computeOptionV== Eigen::ComputeFullV) {
            fullV=true;
        }
    };
    
    SVD_Jb(const Eigen::MatrixXd &A){
        fullU= false;
        fullV= false;
        svd= Eigen::JacobiSVD<Eigen::MatrixXd> (A);
        nrow= A.rows();
        ncol= A.cols();
    };
    virtual ~SVD_Jb(){};
    
    int computeU(){
    // 1 if full U is compted, 0 if thinU is computed, -1 if U is not computed
        if(fullU) return 1;
        else if(svd.computeU()) return 0;
        else return -1;
    };
    int computeV(){
    // 1 if full V is computed, 0 if thin V is computed, -1 if V is not computed
        if(fullV) return 1;
        else if (svd.computeV()) return 0;
        else return -1;
    };
    
    Eigen::MatrixXd matrixU(){
        return svd.matrixU();
    };
    
    Eigen::MatrixXd matrixV(){
        return svd.matrixV();
    };
    Eigen::MatrixXd matrixS(){
        Eigen::VectorXd singularValues= this-> singularValues();
        long size= singularValues.size();
        Eigen::MatrixXd S; S.setIdentity(size,  size);
        S.diagonal()= singularValues;
        
        if(fullU) S.conservativeResize(nrow, size);
        if(fullV) S.conservativeResize(size, ncol);
        
        return S;
    };
    
    Eigen::VectorXd singularValues(){
    // return the singularValues, sorted in decreasing order
        return svd.singularValues();
    };
    
    long nonZeroSingularValues(){
    // return the number of non zero singular values
        return svd.nonzeroSingularValues();
    };
    long rank() {return svd.rank();};
    
    
};

#endif /* SVD_Jb_hpp */
