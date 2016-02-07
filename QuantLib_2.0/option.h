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
    
};

#endif /* option_h */
