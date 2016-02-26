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
#include "options_info.hpp"
#include <Eigen/Dense>
#include <boost/random.hpp>


namespace QLib{
    namespace Simu{
        using namespace QLib;
        using namespace QLib::Simu;
        
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
    SimuEuroBarrier(const barrier_option & o,long path, long time_steps,  unsigned int seed= int(time(0))){
        N= path;
        M= time_steps;
        opt= o;
        
        simuEffective= 1;
        if( opt.barrier>= fmin(opt.K, opt.S) && opt.barrier<= fmax(opt.K, opt.S) ) {
            if (opt.Call && opt.knock_out && opt.barrier<= opt.K) simuEffective= -1;  // knock out call with K>= B >= S. Such an option has no value.
            if (opt.Call && !opt.knock_out && opt.barrier<= opt.K) simuEffective=0; // knock in call with K>=B>= S. Such an option has no difference with a vanilla Euro call.
            if (!opt.Call && !opt.knock_out && opt.barrier>= opt.K ) simuEffective=0; // knock in put with S>=B>=K. Such an option has no difference with a vanilla Euro put.
            if (!opt.Call && opt.knock_out && opt.barrier>= opt.K) simuEffective= -1; // knock out put with S>=B>= K. Such an option has no value.
        }
        
        if (simuEffective==-1) {
            // option value =0
            count_effective= NAN;
            mean= NAN;
            stdiv=NAN;
            N= NAN;
            M= NAN;
            
        } else if(simuEffective==0) {
            // the barrier option can be reduced to vanilla Euro option
            SimuEuro simu_Euro(opt, N, seed);
            M= 1;
            count_effective= N;
            mean= simu_Euro.valuation();
            stdiv= simu_Euro.valuation_stdiv();
            option_value= simu_Euro.optionValueDist();
            asset_path.resize(N, 1);
            asset_path.col(0)= simu_Euro.assetPriceDist();
        } else {
            asset_path.resize(N, M);
            asset_path.setZero();
            option_value.resize(N);
            option_value.setZero();
            
            boost::mt19937 eng(seed);
            boost::normal_distribution<double> normal(0.0, 1.0);
            boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>> rng(eng,normal);
            
            double t= opt.T/ M;
            double c= exp((opt.r- opt.q)*t-.5*opt.sigma*opt.sigma*t);
            std::vector<bool> flag(N);// initialized as false, and if barrier gets hit in path j, flag[j]= true
            for (long i= 0; i< N; i++) flag[i]= false;
            
            for (long i=0; i<N ;i++) {
                // initialize the first col of asset_path
                asset_path(i, 0)= opt.S* c* exp(opt.sigma*sqrt(t)* rng());
                flag[i]= ((opt.S> opt.barrier )^ (asset_path(i,0)>opt.barrier)) || (asset_path(i,0)== opt.barrier);
            }
            
            for (long i= 0; i<N; i++) {
                for (long j=1; j<M; j++) {
                    asset_path(i,j)= asset_path(i,j-1)* c* exp(opt.sigma*sqrt(t)* rng());
                    flag[i]= flag[i] ||(asset_path(i,j)== opt.barrier)|| ((asset_path(i,j-1)> opt.barrier) ^(asset_path(i,j)> opt.barrier));
                }
            }
            
            if (opt.knock_out) {
                // knock out option
                count_effective= N- std::accumulate(flag.begin(), flag.end(), 0);
                for (long i= 0; i<N; i++) {
                    option_value(i)= (!flag[i]) * fmax(0, (-1+ 2* opt.Call) *(asset_path(i,M-1)- opt.K));
                }
                
            } else { // knock in option
                count_effective= std::accumulate(flag.begin(), flag.end(), 0);
                for (long i= 0; i<N; i++) {
                    option_value(i)= (flag[i]) * fmax(0, (-1+ 2* opt.Call) *(asset_path(i,M-1)- opt.K));
                }
            }
            
            mean= option_value.sum()/ option_value.size() * exp(-opt.T*opt.r);
            stdiv= option_value.squaredNorm()/ option_value.size()* exp(-opt.r*opt.T *2);
            stdiv= stdiv- pow(mean,2.0);
            stdiv= sqrt(stdiv/ N);
        }

    };
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)
    SimuEuroBarrier(const barrier_option & o, long path, long time_steps, const std::vector<double> & RN){
        N= path;
        M= time_steps;
        opt= o;
        
        simuEffective= 1;
        if( opt.barrier>= fmin(opt.K, opt.S) && opt.barrier<= fmax(opt.K, opt.S) ) {
            if (opt.Call && opt.knock_out && opt.barrier<= opt.K) simuEffective= -1;  // knock out call with K>= B >= S. Such an option has no value.
            if (opt.Call && !opt.knock_out && opt.barrier<= opt.K) simuEffective=0; // knock in call with K>=B>= S. Such an option has no difference with a vanilla Euro call.
            if (!opt.Call && !opt.knock_out && opt.barrier>= opt.K ) simuEffective=0; // knock in put with S>=B>=K. Such an option has no difference with a vanilla Euro put.
            if (!opt.Call && opt.knock_out && opt.barrier>= opt.K) simuEffective= -1; // knock out put with S>=B>= K. Such an option has no value.
        }
        
        if (simuEffective==-1) {
            // option value =0
            count_effective= NAN;
            mean= NAN;
            stdiv=NAN;
            N= NAN;
            M= NAN;
            
        } else if(simuEffective==0) {
            // the barrier option can be reduced to vanilla Euro option
            SimuEuro simu_Euro(opt, N, RN);
            M= 1;
            count_effective= N;
            mean= simu_Euro.valuation();
            stdiv= simu_Euro.valuation_stdiv();
            option_value= simu_Euro.optionValueDist();
            asset_path.resize(N, 1);
            asset_path.col(0)= simu_Euro.assetPriceDist();
        } else {
            asset_path.resize(N, M);
            asset_path.setZero();
            option_value.resize(N);
            option_value.setZero();
            
            std::vector<const double>::iterator it= RN.begin();
            double t= opt.T/ M;
            double c= exp((opt.r- opt.q)*t-.5*opt.sigma*opt.sigma*t);
            std::vector<bool> flag(N);// initialized as false, and if barrier gets hit in path j, flag[j]= true
            for (long i= 0; i< N; i++) flag[i]= false;
            
            for (long i=0; i<N ;i++) {
                // initialize the first col of asset_path
                asset_path(i, 0)= opt.S* c* exp(opt.sigma*sqrt(t)* (*it));it++;
                flag[i]= ((opt.S> opt.barrier )^ (asset_path(i,0)>opt.barrier)) || (asset_path(i,0)== opt.barrier);
            }
            
            for (long i= 0; i<N; i++) {
                for (long j=1; j<M; j++) {
                    asset_path(i,j)= asset_path(i,j-1)* c* exp(opt.sigma*sqrt(t)* (*it)); it++;
                    flag[i]= flag[i] ||(asset_path(i,j)== opt.barrier)|| ((asset_path(i,j-1)> opt.barrier) ^(asset_path(i,j)> opt.barrier));
                }
            }
            
            if (opt.knock_out) {
                // knock out option
                count_effective= N- std::accumulate(flag.begin(), flag.end(), 0);
                for (long i= 0; i<N; i++) {
                    option_value(i)= (!flag[i]) * fmax(0, (-1+ 2* opt.Call) *(asset_path(i,M-1)- opt.K));
                }
                
            } else { // knock in option
                count_effective= std::accumulate(flag.begin(), flag.end(), 0);
                for (long i= 0; i<N; i++) {
                    option_value(i)= (flag[i]) * fmax(0, (-1+ 2* opt.Call) *(asset_path(i,M-1)- opt.K));
                }
            }
            
            
            mean= option_value.sum()/ option_value.size() * exp(-opt.T*opt.r);
            stdiv= option_value.squaredNorm()/ option_value.size()* exp(-opt.r*opt.T *2);
            stdiv= stdiv- pow(mean,2.0);
            stdiv= sqrt(stdiv/ N);
        }

    };
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
    }
}
 


#endif /* SimuEuroBarrier_hpp */

