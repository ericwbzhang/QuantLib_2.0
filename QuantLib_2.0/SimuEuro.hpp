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
    SimuEuro(const option &o, long path, unsigned int seed=  int(time(0))) ;
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)

    SimuEuro(const option & o, long path, const std::vector<double> & RN);
    // RN is a std vector containing the std normal random number to be used in the simulation. Its length must be no smaller than path.
    virtual ~SimuEuro(){};

    Eigen::VectorXd assetPriceDist(){return asset_price;}
    Eigen::VectorXd  optionValueDist() {return option_value;};
    double valuation() {return mean;};
    double valuation_stdiv() {return stdiv; };

};

#endif /* SimuEuro_hpp */
