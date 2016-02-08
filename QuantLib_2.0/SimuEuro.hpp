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
#include "option.h"
#include <boost/numeric/ublas/vector.hpp>

class SimuEuro{
    // we assume the asset price follows log normal distribution, and the parameters come from the opt object.
protected:
    long N;
    option opt;
    boost::numeric::ublas::vector<double> asset_price;
    boost::numeric::ublas::vector<double> option_value;
    long double mean;
    long double variance;
    
public:
    SimuEuro(){};
    SimuEuro(option o, long path, unsigned int seed=  int(time(0))) ;
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)
    
    SimuEuro(option o, long path, std::vector<double> RN);
    // RN is a std vector containing the std normal random number to be used in the simulation. Its length must be no smaller than path.
    virtual ~SimuEuro(){};
    
    boost::numeric::ublas::vector<double > assetPriceDist(){return asset_price;}
    boost::numeric::ublas::vector<double>  optionValueDist() {return option_value;};
    double valuation() {return mean;};
    double stdiv() {return sqrt(variance); };
    
};

#endif /* SimuEuro_hpp */
