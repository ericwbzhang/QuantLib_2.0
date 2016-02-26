//
//  multiNormalRN.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/10/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef multiNormalRN_hpp
#define multiNormalRN_hpp

#include <stdio.h>
#include <boost/random.hpp>
#include <Eigen/Dense>
namespace QLib{
    using namespace QLib;

class multiNormalRN{
    // multiNormal generates multivariate normal distributed RN given the mean vector and covariance matrix. It utlizes the one dimension RN generator in boost (mt19937 random engine).
protected:
    Eigen::VectorXd mean;
    Eigen::MatrixXd cov;
    Eigen::MatrixXd U;
    // cov= U'U;
    
public:
    multiNormalRN(){};
    multiNormalRN( const Eigen::VectorXd & m, const Eigen::MatrixXd & omega){
        mean= m;
        cov= omega;
        Eigen::LLT<Eigen::MatrixXd> chol(cov);
        U= chol.matrixU();
    };
    
    virtual ~multiNormalRN(){};
    
    Eigen::MatrixXd sample(long n, unsigned int seed= int(time(0))){
        long p= mean.size();
        Eigen::MatrixXd X(n, p);
        
        boost::mt19937 eng(seed);
        boost::normal_distribution<double> normal(0.0, 1.0);
        boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>> rng(eng, normal);
        
        Eigen::VectorXd y(p);
        for (long i= 0; i< n; i++) {
            y.setZero();
            for (long j= 0; j< p ; j++) y(j)= rng();
            
            X.row(i)= y.transpose()*U;
        }
        return X;

    };
    // it returns a sample of size n. Each record has p elements, p= cov.size1();
    
    Eigen::MatrixXd sample(long n, const Eigen::MatrixXd & RN){
        Eigen::MatrixXd X=RN;
        for (long i= 0; i< n; i++) {
            X.row(i)= (U.transpose()* X.row(i).transpose()). transpose();
        }
        return X;

    };
    // RN is a matrix of n*p. It contains iid std normal distributed RNs.
    
    Eigen::VectorXd dist_mean(){return mean;};
    Eigen::MatrixXd dist_cov() {return cov; };
    
};
}

#endif /* multiNormalRN_hpp */


