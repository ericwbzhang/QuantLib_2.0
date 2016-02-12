//
//  SimuEuroBarrier.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/8/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef SimuEuroBarrier_hpp
#define SimuEuroBarrier_hpp

#include <stdio.h>
#include "options.hpp"
#include <Eigen/Dense>

class SimuEuroBarrier{
// Euro Barrier option price strongly depends on the underlying asset price path. Once the barrier gets hit, the option value stays at 0.
// We assume the asset price follows geometrical BM and the parameters come from the option_barrier object.
    
protected:
    long N,M;
    barrier_option opt;
    Eigen::MatrixXd asset_path;
    Eigen::VectorXd option_value;
    int simuEffective;
    // Note some barrier option is not real barrier option-- they either have value 0 for sure, or can be reduced to vanilla Euro options.
    // simuEffective = 1 if we do need to run the simulation to price the barrier option; -1 if we really dont need to run any simulation and know that the value is 0 for sure; 0 if we can reduce the barrier option to vanilla Euro option and use simple simulation to solve it.
    
    
    long count_effective;
    // count how many paths which are not knocked out(for knock out option) or are knocked in (for knock in option).
    long double mean;
    long double stdiv;
    
    
public:
    SimuEuroBarrier(){};
    SimuEuroBarrier(const barrier_option & o,long path, long time_steps,  unsigned int seed= int(time(0)));
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)
    SimuEuroBarrier(const barrier_option & o, long path, long time_steps, const std::vector<double> & RN);
    // RN is a std vector containing the std normal random number to be used in the simulation. Its length must be no smaller than N*M
    
    virtual ~SimuEuroBarrier(){};
    
    Eigen::MatrixXd assetPricePath() {return asset_path;}
    Eigen::VectorXd optionValueDist(){return option_value;}
    double effectRatio() {return count_effective/double(N);};
    // effectRatio returns count_effective/N, the fraction of paths which are not knocked out(for knock out option) or are knocked in (for knock in option). It returns NAN if the option does not need simulation pricers.
    
    int simuEffect() {return simuEffective; }
    double valuation() {return mean;}
    double valuation_stdiv(){return stdiv; }
    
};

 


#endif /* SimuEuroBarrier_hpp */

