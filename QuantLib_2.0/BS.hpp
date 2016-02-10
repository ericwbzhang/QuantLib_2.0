//
//  BS.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef BS_hpp
#define BS_hpp

#include <stdio.h>
#include "options.hpp"
#include <boost/math/distributions/normal.hpp>

struct BS{

    option opt;
    double d1 ,d2, value;
    
    BS(){}
    BS(option o){
        opt=o;
        d1=(log(opt.S)-log(opt.K)+(opt.r-opt.q+0.5*opt.sigma*opt.sigma)*opt.T)/(opt.sigma*sqrt(opt.T));
        d2=d1- opt.sigma*sqrt(opt.T);
        
        boost::math::normal_distribution<double> normal(0,1.0);
        
        if (opt.Call==1) //Call
        {
            value= opt.S* exp(-opt.q *opt.T )* cdf(normal, d1)- opt.K* exp(-opt.r* opt.T)* cdf(normal, d2);
            
        }else{  //Put
            value= opt.K* exp(-opt.r*opt.T)* cdf(normal, -d2)- opt.S* exp(-opt.q* opt.T)* cdf(normal, -d1);
        }
        
    }
    
    
    double price(){
        return value;
    }
    double delta(){
        double result;
        boost::math::normal_distribution<double > normal(0,1.0);
        
        if (opt.Call==1){
            // Call
            result= exp(-opt.q*opt.T) * cdf(normal, d1);
            
        }else {
            // Put
            result= -exp( -opt.q*opt.T)* cdf(normal, -d1);
            
        }
        
        return result;
        
    }
    
    double gamma(){
        boost::math::normal_distribution<double> normal(0,1.0);
        
        return exp(-opt.q* opt.T)* cdf(normal, d1)/ (opt.S* opt.sigma* sqrt(opt.T));
        
    }
    
    double omega(){
        return this->delta()* opt.S/this->price();
    }
    
    double psi(){
        double result;
        boost::math::normal_distribution<double> normal(0,1.0);
        
        if (opt.Call==1){
            //Call
            result= -pdf(normal, d1) * sqrt( opt.T)/ opt.sigma- opt.S* cdf( normal, d1)* exp(-opt.q*opt.T)*opt.T+ opt.K*exp(-opt.r*opt.T)* pdf(normal, d2)* sqrt( opt.T)/ opt.sigma;
        }else{
            //Put
            
            result= -pdf(normal, d1)* sqrt(opt.T)/opt.sigma - opt.S* cdf(normal, d1)*exp(-opt.q*opt.T)* opt.T+ opt.K* exp(-opt.r* opt.T)* pdf(normal, d2)*sqrt(opt.T)/ opt.sigma+ opt.S* exp(-opt.q* opt.T)* opt.T;
        
        }
        
        return result;
    }
    
    double rho(){
        double result;
        boost::math::normal_distribution<double > normal(0,1.0);
        
        if (opt.Call==1) {
            //Call
            result= opt.K* opt.T* exp(-opt.r* opt.T) * cdf(normal, d2);
            
        }else {
            //put
            result= -opt.K * opt.T * exp(-opt.r*opt.T)* cdf(normal, -d2);
            
        }
        return result;
    }
    
    double theta(){
        double result;
        boost::math::normal_distribution<double> normal(0,1.0);
        
        if (opt.Call==1){
            //Call
            result= -exp(-opt.q*opt.T) * opt.S* pdf(normal, d1)* opt.sigma/ (2*sqrt(opt.T)- opt.r* opt.K*exp( -opt.r* opt.T)* cdf(normal, d2)+ opt.q* opt.S*exp( -opt.q*opt.T)* cdf(normal, d1));
        }else{
        //Put
            result= -exp(-opt.q*opt.T) * opt.S* pdf(normal, d1)* opt.sigma/ (2*sqrt(opt.T)+ opt.r* opt.K*exp( -opt.r* opt.T)* cdf(normal, -d2)- opt.q* opt.S*exp( -opt.q*opt.T)* cdf(normal, -d1));
        }
        
        return result;
    }
    
    double vega(){
        boost::math::normal_distribution<double> normal(0,1.0);
        return opt.S* exp( -opt.q* opt.T)* pdf(normal, d1)* sqrt(opt.T);
    }
};

#endif /* BS_hpp */
