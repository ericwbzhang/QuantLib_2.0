//
//  SimuEuro.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/8/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "SimuEuro.hpp"
#include <boost/random.hpp>

SimuEuro::SimuEuro(const option & o, long path, unsigned int seed) {

    opt=o;
    N= path;
    asset_price.resize(N);
    asset_price.setZero();
    option_value= asset_price;
    
    boost::mt19937 eng(seed);
    boost::normal_distribution<double> normal(0.0, 1.0);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double>> rng(eng, normal);
    
    for (long i=0; i< N; i++) {
        asset_price(i)=opt.S* exp((opt.r- opt.q)*opt.T-.5*opt.sigma*opt.sigma*opt.T+ opt.sigma* sqrt(opt.T)* rng());
        
        if(opt.Call) option_value(i)= fmax(asset_price(i)- opt.K,0.0);
        else option_value(i)= fmax(-asset_price(i)+opt.K, 0.0);
    }
    
    mean= option_value.sum()/ option_value.size() * exp(-opt.T*opt.r);
    stdiv= option_value.squaredNorm()/ option_value.size()* exp(-opt.r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);
}

SimuEuro::SimuEuro(const option & o, long path, const std::vector<double> & RN){
    opt=o;
    N= path;
    asset_price.resize(N);
    asset_price.setZero();
    option_value= asset_price;
    
    for (long i=0; i< N; i++) {
        asset_price(i)=opt.S* exp((opt.r- opt.q)*opt.T-.5*opt.sigma*opt.sigma*opt.T+ opt.sigma* sqrt(opt.T)* RN[i]);
        
        if(opt.Call) option_value(i)= fmax(asset_price(i)- opt.K,0.0);
        else option_value(i)= fmax(-asset_price(i)+opt.K, 0.0);
    }
    
    mean= option_value.sum()/ option_value.size() * exp(-opt.T*opt.r);
    stdiv= option_value.squaredNorm()/ option_value.size()* exp(-opt.r*opt.T *2);
    stdiv= stdiv- pow(mean,2.0);
    stdiv= sqrt(stdiv/ N);
    
}
