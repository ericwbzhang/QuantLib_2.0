//
//  SimuEuro.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/8/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef SimuEuro_hpp
#define SimuEuro_hpp

#include <stdio.h>
#include "options_info.hpp"
#include <Eigen/Dense>
#include <boost/random.hpp>
namespace QLib{
    namespace Simu{
        using namespace QLib;
        using namespace QLib::Simu;
class SimuEuro{
    // we assume the asset price follows log normal distribution, and the parameters come from the opt object.
protected:
    long N;
    option opt;
    Eigen::VectorXd asset_price;
    Eigen::VectorXd option_value;
    long double mean;
    long double stdiv;
    // mean is the simulation avg
    // stdiv is the simulation avg's stdiv


public:
    SimuEuro(){};
    SimuEuro(const option &o, long path, unsigned int seed=  int(time(0))) {
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

    };
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)

    SimuEuro(const option & o, long path, const std::vector<double> & RN){
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
    };
    // RN is a std vector containing the std normal random number to be used in the simulation. Its length must be no smaller than path.
    virtual ~SimuEuro(){};

    Eigen::VectorXd assetPriceDist(){return asset_price;}
    Eigen::VectorXd  optionValueDist() {return option_value;};
    double valuation() {return mean;};
    double valuation_stdiv() {return stdiv; };

};
    }
}

#endif /* SimuEuro_hpp */
