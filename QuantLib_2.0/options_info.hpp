//
//  options_info.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef options_hpp
#define options_hpp

#include <stdio.h>
#include <math.h>
#include <vector>
#include <Math/IFunction.h>
#include <Eigen/Dense>

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

    barrier_option(const option & o, double bar, int k_out): option(o), barrier(bar), knock_out(k_out){};

    virtual ~barrier_option(){};

};

struct asset{

    //asset contains the parameters of a price process. We assume the price follows log normal distribution.

    double S, r, q, sigma;

    asset(){};
    asset(double stock_price, double rf, double div, double volatility) : S(stock_price), r(rf), q(div), sigma(volatility){};
    asset(const option & o) : S(o.S), r(o.r), q(o.q), sigma(o.sigma){};

    virtual ~asset(){};

    double expected_price(double T) {return S* exp((r-q)*T);};

};

struct nonPathDependentBasket_option{
    // nonPathDependentBasket_option contains the option parameters for a non path dependent basket option.
    // The payoff is only related to the asset price at final stage. It takes form as a real value function from R^p to R, where p is the number of assets we are considering.
    // We assume the assetp price follows log normal distribution, the distribution parameters are specified by the option params.
    std::vector<asset> asset_vec;
    long count_assets;
    double T;
    std::shared_ptr<ROOT::Math::IBaseFunctionMultiDim> f;
    Eigen::MatrixXd cov;


    nonPathDependentBasket_option(){};
    nonPathDependentBasket_option(const std::vector<asset> & vec, double expiration, ROOT::Math::IBaseFunctionMultiDim * payoff, const Eigen::MatrixXd & mtx): asset_vec(vec), T(expiration), cov(mtx), count_assets(vec.size()){
        f.reset( payoff->Clone());
    };
    nonPathDependentBasket_option( const nonPathDependentBasket_option & o): asset_vec(o.asset_vec), count_assets(o.count_assets), T(o.T), cov(o.cov){
        // copy constructor. we want a deep copy, ie, we also need to clone *f

        f.reset(o.f->Clone());

    };

    virtual ~nonPathDependentBasket_option () {};

    nonPathDependentBasket_option& operator =( const nonPathDependentBasket_option & o) {

        asset_vec= o.asset_vec;
        count_assets= o.count_assets;
        T=o.T;
        cov=o.cov;

        f.reset(o.f->Clone());

        return *this;
    };

    std::vector<double > expected_price(){
        std::vector<double> p(count_assets);
        for(long i=0; i< count_assets; i++) {
            p[i]= asset_vec[i].expected_price(T);
        }

        return p;
    };

    Eigen::MatrixXd corr(){

        long nrow= cov.rows();
        long ncol= cov.cols();
        Eigen::MatrixXd Rho(nrow, ncol); Rho.setZero();


        for (long i= 0; i< nrow; i++)
            for ( long j=0; j< ncol; j++)
                Rho(i,j)= cov(i,j)/ sqrt(cov(i,i)* cov(j,j));

        return Rho;
    };
};



#endif /* options_hpp */
