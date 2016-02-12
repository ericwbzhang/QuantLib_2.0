//
//  BinomialBSTree.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "BinomialBSTree.hpp"

BinomialBSTree:: BinomialBSTree(const option &o, long steps){
    N= steps;
    opt= o;
    
    t= opt.T/ N;
    u= exp(opt.sigma* sqrt(t));
    d= 1/u;
    p= (exp((opt.r-opt.q)*t)-d)/(u-d);
    q=1-p;
    disc= exp(-opt.r*t);
    
    value.resize(N);
    value.setZero();
    
    // Evolve the value vector
    
    option o1= opt;
    o1.T= t;
    if (opt.Euro) {
        // For Euro option, the call and put share same evolution algorithm.
        for(long i=0; i< N; i++) {
            o1.S= opt.S* pow(u, N-1-i)* pow(d, i);
            BS bs(o1);
            value(i)= bs.price();
        }
        
        for (long j=0; j<N-1; j++) {
            for (long i=0; i<N-1-j ;i++) {
                value(i)= disc* (value(i)*p + value(i+1)* q);
            }
        }
        
    }else {
        if (opt.Call) {// American Call
            for (long i=0; i< N; i++) {
                o1.S= opt.S* pow(u, N-1-i)* pow(d, i);
                BS bs(o1);
                value(i)= bs.price();
                value(i)= fmax( value(i), o1.S- opt.K);
            }
            
            for (long j=0; j< N-1; j++) {
                for(long i=0; i< N-1-j; i++){
                    value(i)= disc* (p*value(i)+ q* value(i+1));
                    value(i)= fmax(value(i), opt.S* pow(d, i)* pow(u, N-2-j-i)- opt.K);
                    
                }
            }
        }else { //American Put
            for (long i=0; i<N; i++) {
                o1.S= opt.S* pow(u, N-1-i)* pow(d, i);
                BS bs(o1);
                value(i)= bs.price();
                value(i)= fmax( value(i), opt.K-o1.S);
            }
            
            for (long j=0; j< N-1; j++) {
                for(long i=0; i< N-1-j; i++){
                    value(i)= disc* (p*value(i)+ q* value(i+1));
                    value(i)= fmax(value(i), -opt.S* pow(d, i)* pow(u, N-2-j-i)+ opt.K);
                    
                }
            }
        }
    }
    
}
