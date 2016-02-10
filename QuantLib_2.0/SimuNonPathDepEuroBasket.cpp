//
//  SimuNonPathDepEuroBasket.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "SimuNonPathDepEuroBasket.hpp"
#include <boost/random.hpp>
#include "multiNormal.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>


SimuNonPathDepEuroBasket::SimuNonPathDepEuroBasket(nonPathDependentBasket_option o, long paths, unsigned int seed){
    
    N=paths;
    opt=o;
    option_value.resize(N); option_value.clear();
    asset_price.resize(N, opt.count_assets); asset_price.clear();
    
    // Boost has no RN generator for multivariate normal. I implement it, in class multiNormal. It involves Cholesky decomp to transform a vector of iid std normal RN to a sample of multivariate normal distribution with given covariance matrix.
    boost::numeric::ublas::vector<double> m(opt.count_assets); m.clear();
    multiNormal multiNormalRNG(m, opt.cov);
    asset_price= multiNormalRNG.sample(N, seed);
    
    for (long j=0; j<opt.count_assets; j++) {
        asset a= opt.asset_vec[j];
        long double c= a.expected_price(opt.T) * exp(-.5* a.sigma*a.sigma* opt.T);
        for (long i=0; i<N; i++){
            asset_price(i,j)= c* exp(sqrt(opt.T)* asset_price(i,j));
        }
    }
    
    for (long i=0; i<N; i++){
        option_value(i)= opt.f->operator()(boost::numeric::ublas::row(asset_price, i));
    }
    
    
    double r= opt.asset_vec[0].r;
    
    mean= boost::numeric::ublas::sum(option_value)/ option_value.size() * exp(-opt.T*r);
    stdiv= pow(boost::numeric::ublas::norm_2(option_value), 2.0)/ option_value.size()* exp(-r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);
}


SimuNonPathDepEuroBasket::SimuNonPathDepEuroBasket(nonPathDependentBasket_option o, long paths, boost::numeric::ublas::matrix<double> RN){
    
    N=paths;
    opt=o;
    option_value.resize(N); option_value.clear();
    asset_price.resize(N, opt.count_assets); asset_price.clear();
    
    // Boost has no RN generator for multivariate normal. I implement it, in class multiNormal. It involves Cholesky decomp to transform a vector of iid std normal RN to a sample of multivariate normal distribution with given covariance matrix.
    boost::numeric::ublas::vector<double> m(opt.count_assets); m.clear();
    multiNormal multiNormalRNG(m, opt.cov);
    asset_price= multiNormalRNG.sample(N, RN);
    
    for (long j=0; j<opt.count_assets; j++) {
        asset a= opt.asset_vec[j];
        long double c= a.expected_price(opt.T) * exp(-.5* a.sigma*a.sigma* opt.T);
        for (long i=0; i<N; i++){
            asset_price(i,j)= c* exp(sqrt(opt.T)* asset_price(i,j));
        }
    }
    
    for (long i=0; i<N; i++){
        option_value(i)= opt.f->operator()(boost::numeric::ublas::row(asset_price, i));
    }
    
    
    double r= opt.asset_vec[0].r;
    
    mean= boost::numeric::ublas::sum(option_value)/ option_value.size() * exp(-opt.T*r);
    stdiv= pow(boost::numeric::ublas::norm_2(option_value), 2.0)/ option_value.size()* exp(-r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);

    
}