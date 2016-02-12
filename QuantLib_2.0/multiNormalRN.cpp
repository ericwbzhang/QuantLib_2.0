//
//  multiNormalRN.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/10/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "multiNormalRN.hpp"
#include "algebra.h"

multiNormalRN::multiNormalRN(const Eigen::VectorXd &m, const Eigen::MatrixXd &omega) {
    mean= m;
    cov= omega;
    Cholesky_LLT chol(cov);
    U= chol.factorU();
}

Eigen::MatrixXd multiNormalRN::sample(long n, unsigned int seed) {
    
    long p= mean.size();
    Eigen::MatrixXd X(n, p);
    
    boost::mt19937 eng(seed);
    boost::normal_distribution<double> normal(0.0, 1.0);
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>> rng(eng, normal);
    
    Eigen::VectorXd y(p);
    for (long i= 0; i< n; i++) {
        y.setZero();
        for (long j= 0; j< p ; j++) y(j)= rng();
        
        X.row(i)= (U.transpose()* y).transpose();
    }
    return X;
}


Eigen::MatrixXd multiNormalRN::sample(long n, const Eigen::MatrixXd & RN){
    
    Eigen::MatrixXd X=RN;
    for (long i= 0; i< n; i++) {
        X.row(i)= (U.transpose()* X.row(i).transpose()). transpose();
    }
    return X;
    
}
