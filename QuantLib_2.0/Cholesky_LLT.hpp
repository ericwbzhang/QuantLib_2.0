//
//  Cholesky_LLT.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/11/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef Cholesky_LLT_hpp
#define Cholesky_LLT_hpp


#include <stdio.h> 
#include <Eigen/Dense> 

class Cholesky_LLT {

protected:
    Eigen::MatrixXd U;

public:
    Cholesky_LLT (){};
    Cholesky_LLT (const Eigen::MatrixXd &A){
        U= A.llt().matrixU();
    }
    
    virtual ~Cholesky_LLT(){};
    
    Eigen::MatrixXd factorU(){return U;};
    
};

#endif /* Cholesky_LLT_hpp */
