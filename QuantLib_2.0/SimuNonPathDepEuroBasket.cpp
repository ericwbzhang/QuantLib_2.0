//
//  SimuNonPathDepEuroBasket.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "SimuNonPathDepEuroBasket.hpp"
#include <boost/random.hpp>
#include "multiNormalRN.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>


SimuNonPathDepEuroBasket::SimuNonPathDepEuroBasket(const nonPathDependentBasket_option & o, long paths, unsigned int seed){
    
    N=paths;
    opt=o;
    option_value.resize(N); option_value.setZero();
    asset_price.resize(N, opt.count_assets); asset_price.setZero();
    
    // Boost has no RN generator for multivariate normal. I implement it, in class multiNormal. It involves Cholesky decomp to transform a vector of iid std normal RN to a sample of multivariate normal distribution with given covariance matrix.
    Eigen::VectorXd m(opt.count_assets); m.setZero();
    multiNormalRN multiNormalRNG(m, opt.cov);
    asset_price= multiNormalRNG.sample(N, seed);
    
    for (long j=0; j<opt.count_assets; j++) {
        asset a= opt.asset_vec[j];
        long double c= a.expected_price(opt.T) * exp(-.5* a.sigma*a.sigma* opt.T);
        for (long i=0; i<N; i++){
            asset_price(i,j)= c* exp(sqrt(opt.T)* asset_price(i,j));
        }
    }
    
    for (long i=0; i<N; i++){
//        Eigen::VectorXd t(asset_price.row(i));
//        std::vector<double> tmp(opt.count_assets);
//        for (long i=0; i< tmp.size(); i++) tmp[i]= t(i);
        
        option_value(i)= opt.f->operator()(asset_price.row(i));
    }
    
    
    double r= opt.asset_vec[0].r;
    
    
    mean= option_value.sum()/ option_value.size() * exp(-opt.T*r);
    stdiv= option_value.squaredNorm()/ option_value.size()* exp(-r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);
}


SimuNonPathDepEuroBasket::SimuNonPathDepEuroBasket(const nonPathDependentBasket_option & o, long paths, const Eigen::MatrixXd & RN){
    
    N=paths;
    opt=o;
    option_value.resize(N); option_value.setZero();
    asset_price.resize(N, opt.count_assets); asset_price.setZero();
    
    // Boost has no RN generator for multivariate normal. I implement it, in class multiNormal. It involves Cholesky decomp to transform a vector of iid std normal RN to a sample of multivariate normal distribution with given covariance matrix.
    Eigen::VectorXd m(opt.count_assets); m.setZero();
    multiNormalRN multiNormalRNG(m, opt.cov);
    asset_price= multiNormalRNG.sample(N, RN);
    
    for (long j=0; j<opt.count_assets; j++) {
        asset a= opt.asset_vec[j];
        long double c= a.expected_price(opt.T) * exp(-.5* a.sigma*a.sigma* opt.T);
        for (long i=0; i<N; i++){
            asset_price(i,j)= c* exp(sqrt(opt.T)* asset_price(i,j));
        }
    }
    
    for (long i=0; i<N; i++){
//        Eigen::VectorXd t(asset_price.row(i));
//        std::vector<double> tmp(opt.count_assets);
//        for (long i=0; i< tmp.size(); i++) tmp[i]= t(i);
        
        option_value(i)= opt.f->operator()(asset_price.row(i));
    }
    
    
    double r= opt.asset_vec[0].r;
    
    mean= option_value.sum()/ option_value.size() * exp(-opt.T*r);
    stdiv= option_value.squaredNorm()/ option_value.size()* exp(-r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);
    
}