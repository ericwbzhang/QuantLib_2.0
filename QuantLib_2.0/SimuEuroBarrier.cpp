//
//  SimuEuroBarrier.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/8/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "SimuEuroBarrier.hpp"
#include <boost/random.hpp>
#include "SimuEuro.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>



SimuEuroBarrier::SimuEuroBarrier(barrier_option o, long paths, long time_steps, std::vector<double> RN){
    
    N= paths;
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
        boost::numeric::ublas::column(asset_path, 0)= simu_Euro.assetPriceDist();
    } else {
        asset_path.resize(N, M);
        asset_path.clear();
        option_value.resize(N);
        option_value.clear();
        
        std::vector<double>::iterator it= RN.begin();
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
        
        mean= boost::numeric::ublas::sum(option_value)/ option_value.size() * exp(-opt.T*opt.r);
        stdiv= pow(boost::numeric::ublas::norm_2(option_value), 2.0)/ option_value.size()* exp(-opt.r*opt.T *2);
        stdiv= stdiv- pow(mean,2.0);
        stdiv= sqrt(stdiv/ N);
    }
    
}


SimuEuroBarrier::SimuEuroBarrier(barrier_option o, long paths, long time_steps, unsigned int seed){
    
    N= paths;
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
        boost::numeric::ublas::column(asset_path, 0)= simu_Euro.assetPriceDist();
    } else {
        asset_path.resize(N, M);
        asset_path.clear();
        option_value.resize(N);
        option_value.clear();
        
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
        
        mean= boost::numeric::ublas::sum(option_value)/ option_value.size() * exp(-opt.T*opt.r);
        stdiv= pow(boost::numeric::ublas::norm_2(option_value), 2.0)/ option_value.size()* exp(-opt.r*opt.T *2);
        stdiv= stdiv- pow(mean,2.0);
        stdiv= sqrt(stdiv/ N);
    }


}




