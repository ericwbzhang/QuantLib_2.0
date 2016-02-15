//
//  SimuLS.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/14/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "SimuLS.hpp"
#include <boost/random.hpp>

SimuLS::SimuLS(const option& o, long path, long time_steps, long effpath_lowLimit, unsigned int seed){

    opt=o;
    N= path;
    M= time_steps;
    effpathcountMIN= effpath_lowLimit;
    
    if( opt.Euro) { // Euro option
        SimuEuro SE(opt,N);
        asset_path= SE.assetPriceDist();
        option_cashflow= SE.optionValueDist();
        
        simuEffective= 0;
        earlyExer_count=NAN;
        mean= SE.valuation();
        stdiv= SE.valuation_stdiv();
        
    }else{ //American option
        simuEffective=1;
        asset_path.setZero(N, M);
        option_cashflow. setZero(N, M);
        //  evaluate asset_path
        boost::mt19937 eng(seed);
        boost::normal_distribution<double> normal(0.0, 1.0);
        boost::variate_generator<boost::mt19937 & , boost::normal_distribution<>> rng(eng, normal);
        double t= opt.T/ time_steps;
        double c= exp((opt.r- opt.q)*t-.5*opt.sigma*opt.sigma*t);
        double sigma_sqrtt= opt.sigma* sqrt(t);
        
        
        for (long i=0; i<N; i++) {
            asset_path(i,0)= opt.S* c* exp(sigma_sqrtt* rng()) ;
        }
        
        for (long i= 0; i<N; i++) {
            for (long j=1; j<M; j++) {
                asset_path(i,j)= asset_path(i,j-1)* c* exp(sigma_sqrtt* rng());
            }
        }
        
        this->compute();
    }
}

SimuLS::SimuLS(const option &o, long path, long time_steps,  const std::vector<double> &RN, long effpath_lowLimit){
    
    opt=o;
    N= path;
    M= time_steps+1;
    effpathcountMIN= effpath_lowLimit;
    
    if( opt.Euro) { // Euro option
        SimuEuro SE(opt,N);
        asset_path= SE.assetPriceDist();
        option_cashflow= SE.optionValueDist();
        
        simuEffective= 0;
        earlyExer_count=NAN;
        mean= SE.valuation();
        stdiv= SE.valuation_stdiv();
        
    }else{ //American option
        simuEffective=1;
        asset_path.setZero(N, M);
        option_cashflow. setZero(N, M);
        //  evaluate asset_path

        
        double t= opt.T/ time_steps;
        double c= exp((opt.r- opt.q)*t-.5*opt.sigma*opt.sigma*t);
        double sigma_sqrtt= opt.sigma* sqrt(t);
        
        long q=0;
        for (long i=0;i <N; i++) {
            asset_path(i,0)= opt.S* c* exp(sigma_sqrtt*RN[q]); q++;
        }
        for (long i= 0; i<N; i++) {
            for (long j=1; j<M; j++) {
                asset_path(i,j)= asset_path(i,j-1)* c* exp(sigma_sqrtt* RN[q]); q++;
            }
        }
        
        this->compute();
    }
}


void SimuLS::compute(){
    
    double t= opt.T/ M;
    double disc= exp(-opt.r* t);
    // evaluate option cash flow
    for (long i=0; i<N; i++) {
        option_cashflow(i, M-1)= opt.Call ? fmax(asset_path(i,M-1)- opt.K, 0.0) : fmax(opt.K- asset_path(i,M-1), 0.0);
    }
    
    for(long j=M-2; j>-1; j--) {
        Eigen::MatrixXd asset_price;
        Eigen::MatrixXd disc_optionCF;
        Eigen::MatrixXd ExerValue;
        std::unordered_map<long , long > map;
        this-> effpath(asset_path.col(j), option_cashflow.block(0, j+1, N, M-j-1), disc, asset_price, disc_optionCF, ExerValue,  map);
        // we have  asset_price disc_optionCF as vector/matrix and we have expanded asset_price to second order
        // The next step is to construct the regression formula.
        
        if(asset_price.rows()){
            // We do OLS regression
            if (ExerValue.size()<= effpathcountMIN) simuEffective= -1;
            else{
                lm_OLS lm(asset_price, disc_optionCF);
                // we get the fitted value and compare it with the Exercise return
                Eigen::MatrixXd diff(ExerValue);
                diff-= lm.prediction(asset_price);
                for (long l= 0; l< diff.size();l++){
                    if (diff(l)>0 ) {
                        // exercise option
                        option_cashflow.block(map[l], j+1, 1, M-j-1).setZero();
                        option_cashflow(map[l], j)=ExerValue(l);
                    }
                }
            
            }
        }
    }
    
    earlyExer_count=0;
    for (long i=0; i< option_cashflow.rows(); i++){
        for (long j=0; j< option_cashflow.cols()-1; j++){
            if (option_cashflow(i,j)!=0 ) earlyExer_count++;
        }
    }
    
    Eigen::VectorXd DiscVec (option_cashflow.cols());
    for(long i=0; i< DiscVec.size(); i++) {
        DiscVec(i)= pow(disc, i+1);
    }
    optionValue_dist= option_cashflow* DiscVec;
    mean= optionValue_dist.mean();
    stdiv=  optionValue_dist.squaredNorm()- mean*mean* optionValue_dist.size();
    stdiv= sqrt(stdiv/ (optionValue_dist.size()-1))/ sqrt(N);


}

void SimuLS::effpath(const Eigen::VectorXd &S, const Eigen::MatrixXd & optionCF, double disc, Eigen::MatrixXd & asset_price, Eigen::MatrixXd & disc_optionCF, Eigen::MatrixXd & ExerValue, std::unordered_map<long , long > & map){
    std::vector<double> assetPrice, disc_optionCashFlow, exerValue;
    
    // initialize the disc vector
    long k= optionCF.cols();
    Eigen::VectorXd disc_vec(k);
    for(long j= 0; j<k ; j++) disc_vec(j)= pow(disc, j+1);
    
    // compute  assetPrice and disc_optionCashFlow
    long n= S.rows();
    long m=0;
    for (long i=0; i<n; i++) {
        if ((opt.Call && S(i)> opt.K)|| (! opt.Call && S(i)< opt.K) ) {
            //  Exercise is possible
            assetPrice.push_back(S(i));
            exerValue.push_back(std::abs(opt.K- S(i)));  // note: abs return int, here we should use std::abs which returns double
            disc_optionCashFlow.push_back(optionCF.row(i).dot(disc_vec));
            map.insert({m,i});
            m++;
        }
    }
    
    // Map assetPrice and disc_optionCashFlow to Matrix type
    asset_price.resize(assetPrice.size(), 3);
    asset_price.col(0).setOnes();
    asset_price.col(1)= Eigen::Map<Eigen::MatrixXd> (assetPrice.data(), assetPrice.size(), 1);
    asset_price.col(2)= asset_price.col(1).asDiagonal()* asset_price.col(1);
    
    ExerValue= Eigen::Map<Eigen::MatrixXd> (exerValue.data(), exerValue.size(), 1);
    disc_optionCF= Eigen::Map<Eigen::MatrixXd> (disc_optionCashFlow.data(), disc_optionCashFlow.size(),1);

}