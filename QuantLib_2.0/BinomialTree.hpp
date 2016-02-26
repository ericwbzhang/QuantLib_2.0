//
//  BinomialTree.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef BinomialTree_hpp
#define BinomialTree_hpp

#include <stdio.h>
#include "options_info.hpp"
#include <Eigen/Dense>

namespace QLib{
    namespace Tree{
        using namespace QLib;
        using namespace QLib::Tree;
        
class BinomialTree{
protected:
    long N; // number of time steps
    option opt; // the option to be priced
    Eigen::VectorXd value; // the boost vector to capture the option value evolution
    double t, p, q, u, d, disc;
    // t= length of each time step;
    // p= the risk neutral prob of up; q=1-p, the risk neutral prob of down
    // u= the risk neutral multiplier of up; d the risk neutral multiplier of down
    // disc= the risk neutral discounting factor for one time step

public:
    BinomialTree(){};
    BinomialTree(const option & o, long steps){
        opt=o;
        N= steps;
        
        t= opt.T/ N;
        u= exp(opt.sigma* sqrt(t));
        d= 1/u;
        p= (exp((opt.r-opt.q)*t)-d)/(u-d);
        q=1-p;
        disc= exp(-opt.r*t);
        
        value.resize(N+1);
        value.setZero();
        
        // Evolve the value vector
        if (opt.Euro==1){
            if (opt.Call==1){ // Euro Call
                
                for (long i=0; i<N+1; i++){
                    value(i)= fmax(opt.S* pow(d, i)* pow(u, N-i)- opt.K, 0.0);
                }
                for (long j=0; j<N; j++){
                    for (long i= 0; i<N-j; i++){
                        value(i)= disc* (value(i)*p + value(i+1)*q);
                    }
                }
            }else { // Euro Put
                
                for (long i=0; i<N+1; i++){
                    value(i)= fmax(opt.K- opt.S* pow(d,i)* pow(u,N-i), 0.0);
                }
                for (long j=0; j<N; j++){
                    for (long i= 0; i<N-j; i++){
                        value(i)= disc* (value(i)*p + value(i+1)*q);
                    }
                }
            }
            
        }else {
            if( opt.Call==1){// American Call
                
                for (long i=0; i< N+1; i++){
                    value(i)= fmax(opt.S* pow(d,i)* pow(u, N-i)- opt.K, 0.0);
                }
                for (long j=0; j<N; j++){
                    for (long i=0; i< N-j; i++){
                        value(i)= disc*(p*value(i)+ q* value(i+1));
                        value(i)= fmax( value(i), opt.S* pow( d, i)* pow(u, N-j-i-1)-opt.K);
                    }
                }
            }else {// American Put
                
                for (long i=0; i< N+1; i++){
                    value(i)= fmax(-opt.S* pow(d,i)* pow(u, N-i)+ opt.K, 0.0);
                }
                for (long j=0; j<N; j++){
                    for (long i=0; i< N-j; i++){
                        value(i)= disc*(p*value(i)+ q* value(i+1));
                        value(i)= fmax( value(i), -opt.S* pow( d, i)* pow(u, N-j-i-1)+ opt.K);
                    }
                }
            }
        }
        
    };
    

    
    virtual ~BinomialTree(){};
    
    virtual double valuation(){return value(0);};
    virtual double delta(){
        
        double v1= value(1);
        double v0= (value(0)/disc- q*value(1))/p;
        
        return (v0-v1)/ (opt.S*(u-d));
    };
    virtual double gamma(){double v2= value(2);
        double v1= (value(1)/disc - q* value(2))/p;
        double v0= (value(0)/disc- q*value(1))/p;
        v0= (v0/disc - q*v1)/p;
        
        double a ,b,c, result;
        
        a= (v0-v1)/ (opt.S*(u*u- u*d));
        b= (v1- v2)/ (opt.S*(u*d- d*d));
        
        c= .5* (opt.S*(u*u - u*d));
        
        result= (a-b)/c;
        return result;
    };
    virtual double theta(){
        double v1= (value(1)/disc - q* value(2))/p;
        
        return (v1-value(0))/ (2*t);
        

    };

};
    }
}
#endif /* BinomialTree_hpp */
