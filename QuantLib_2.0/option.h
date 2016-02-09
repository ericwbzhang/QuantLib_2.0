//
//  option.h
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef option_h
#define option_h


struct option{
    double S, K, T, r, q, sigma;
    int Call;// Call option if Call=1; Put if Call=0
    int Euro; // Euro optoin if Euro =1; American if Euro =0;
    
    option(){};
    option(double stock_price, double strike, double expiration, double rf, double div, double volatility, int CorP, int Euro_Amer ): S(stock_price), K(strike), T(expiration), r(rf), q(div), sigma(volatility), Call(CorP), Euro(Euro_Amer) {};
    
    virtual ~option(){};
    
};

struct barrier_option : public option{
    // barrier_option derives from option class, and contains the option parameters including barrier.
    // we consider the basic barrier options, knock_out and knock_in. And we assume that all the barrier options comes with the chance to be activated, ie for example no knock out option with barrier has been knocked.
    
    double barrier;
    int knock_out; // knock_out=1 if the option is a knock_out; knock_out =0 if it is knock_in;
    // we define knock out as once we touch the barrier the option get void; knock in as the once we touch the barrier the option gets activated from void.
    
    barrier_option(){};
    barrier_option(double stock_price, double strike, double expiration, double rf, double div, double volatility, int CorP, int Euro_Amer, double bar, int k_out):option(stock_price, strike, expiration, rf, div, volatility, CorP, Euro_Amer), barrier(bar), knock_out(k_out){};
    
    barrier_option(option o, double bar, int k_out): option(o), barrier(bar), knock_out(k_out){};
    
    virtual ~barrier_option(){};
    
};




#endif /* option_h */
