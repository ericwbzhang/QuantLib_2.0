//
//  SimuLS.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/14/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef SimuLS_hpp
#define SimuLS_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include "options.hpp"
#include "SimuEuro.hpp"
#include "lm.hpp"
#include <unordered_map>



class SimuLS{
// SimuLS implements Least Square MC simulation method for American option pricing.
    
protected:
    long N, M;
    option opt;
    Eigen::MatrixXd asset_path;
    Eigen::MatrixXd option_cashflow;
    Eigen::VectorXd optionValue_dist;
    
    long effpathcountMIN;
    int simuEffective;
    // Note: Simulation method does not do well for deep out of the money optoin since few paths knock the strike. We set a effpathcountMIN as the lower threshold for realiable result. We do regression only when the sample size is greater than effpathcountMIN.
    // simuEffective = 1 if we are pricing American options and effpathcountMIN is never hit; 0 if we are pricing Euro options; -1 if are pricing American options and effpathcountMIN is hit.
    
    long earlyExer_count;
    // earlyExer count how many paths where early exercise happens.
    long double mean;
    long double stdiv;
    
    virtual void effpath(const Eigen::VectorXd &S, const Eigen::MatrixXd & optionCF, double disc, Eigen::MatrixXd & asset_price, Eigen::MatrixXd & disc_optionCF, Eigen::MatrixXd & ExerValue, std::unordered_map<long , long > & map);
    
    virtual void compute();
    
public:
    SimuLS(){};
    virtual ~SimuLS(){};
    SimuLS(const option & o,long path, long time_steps,  long effpath_lowLimit=30, unsigned int seed= int(time(0)));
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)
    SimuLS(const option & o, long path, long time_steps,  const std::vector<double> & RN,long effpath_lowLimit= 30);
    // RN is a std vector containing the std normal random number to be used in the simulation. Its length must be no smaller than N*M
    
    virtual Eigen::MatrixXd assetPricePath() {return asset_path;};
    virtual Eigen::MatrixXd optionCashFlow(){return option_cashflow;};
    virtual Eigen::VectorXd optionValueDist(){return optionValue_dist;};
    virtual double earlyExerRatio() {return earlyExer_count/double(N);};
    // earlyExerRatio returns earlyExer_count/N, the fraction of paths where early exercise happens. It returns NAN if the option is European style.
    
    virtual int simuEffect() {return simuEffective; };
    virtual double valuation() {return mean;};
    virtual double valuation_stdiv(){return stdiv; };
    
};


class sampleCalculator {

    // sampleCalculator is a base class to generate derivatives which can take simulation path matrix as input, and return simulation method pricing.
    
protected:
    double theoraticalValue;
public:
    sampleCalculator(){};
    sampleCalculator(double p) :theoraticalValue(p){};
    virtual ~sampleCalculator(){};
    
    virtual double theoraticalvalue() const  {return theoraticalValue; } ;
    virtual Eigen::VectorXd  valuation_dist(const Eigen::MatrixXd & simulationPath) const =0;
};

class sampleCalculator_Euro: public sampleCalculator{
// sampleCalculator_Euro implements the case for European options
    
protected:
    option opt;
    
public:
    
    sampleCalculator_Euro(){};
    virtual ~sampleCalculator_Euro(){};
    sampleCalculator_Euro( const option& o, double p): sampleCalculator(p) {opt=o;  };
    
    virtual Eigen::VectorXd valuation_dist(const Eigen::MatrixXd & simulationPath) const {
        Eigen::VectorXd res= simulationPath.col(simulationPath.cols()-1);
        for(long i=0; i< res.size(); i++){
            res(i)= opt.Call? fmax(res(i)- opt.K, 0.0) : fmax(opt.K- res(i), 0.0);
        }
        return res* exp(-opt.r* opt.T);
    };
    
};

class SimuLS_CV{
// SimuLS_CV applies control variate method to SimuLS
    
protected:
    
    SimuLS simls;
    double calibration_coef;
    double rho;
    double CV_TheoraticalValue;
    double CV_SampleValue;
    //Eigen::VectorXd optValueDist_calibrated;
    
    void compute(const sampleCalculator & SC){
        Eigen::VectorXd t1= simls.optionValueDist();
        Eigen::VectorXd t2= SC.valuation_dist(simls.assetPricePath());
        
        double m1, m2, var1, var2;
        m1= t1.mean();
        m2= t2.mean();
        var2= t2.squaredNorm()- m2*m2* t2.size(); var2/= (t2.size()-1);
        var1= t1.squaredNorm()- m1*m1* t1.size(); var1/= (t1.size()-1);
        
        calibration_coef= t1.dot(t2)/ t1.size() - m1*m2;
        rho= calibration_coef/sqrt(var1*var2);
        calibration_coef/= -var2;
        
        CV_TheoraticalValue= SC.theoraticalvalue();
        CV_SampleValue= m2;

    };
    
public:
    SimuLS_CV(){};
    virtual ~SimuLS_CV(){};
    SimuLS_CV( const option& o, long path, long time_steps,  const sampleCalculator &  SC, long effpath_lowLimit=30, unsigned int seed= int(time(0))){
        simls= SimuLS(o, path, time_steps, effpath_lowLimit, seed);
        
        this-> compute(SC);
    };
    SimuLS_CV(const option &o,  long path, long time_steps,   const sampleCalculator & SC,  const std::vector<double> & RN,long effpath_lowLimit=30){
        simls= SimuLS(o, path, time_steps, RN, effpath_lowLimit);
        
        this->compute(SC);
    };
    
    
    virtual Eigen::MatrixXd assetPricePath() {return simls.assetPricePath();};
    virtual Eigen::MatrixXd optionCashFlow(){return simls.optionCashFlow();};
    
    
    virtual Eigen::VectorXd optionValueDist_raw(){return simls.optionValueDist();};
    virtual double calibrationCoef() {return calibration_coef; };
    virtual double CVcorr(){return rho; };
    virtual double earlyExerRatio() {return simls.earlyExerRatio();};
    // earlyExerRatio returns earlyExer_count/N, the fraction of paths where early exercise happens. It returns NAN if the option is European style.
    
    virtual int simuEffect() {return simls.simuEffect(); };
    virtual double valuation_raw() {
        return simls.valuation();
    };
    
    virtual double valuation_calibrated(){
        return this->valuation_raw()+ calibration_coef* ( CV_SampleValue- CV_TheoraticalValue);
    };
    virtual double valuation_stdiv_raw(){
        return simls.valuation_stdiv();
    };
    
    virtual double valuation_stdiv_calibrated(){
        return this->valuation_stdiv_raw()* sqrt(1- rho*rho);
    
    };
    
    virtual double stdiv_shrinkage(){return sqrt(1-rho*rho);};
};


/*
 ************ test ************

int main(){
    
    option opt(36,40, 1, 0.06, 0, 0.4, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
    
    opt.Call =0;
    BS bs(opt);
    std:: cout<< bs.price()<< std::endl;
    
    opt.Euro=0;
    
    boost::mt19937 eng(time(0));
    boost::normal_distribution<> normal(0,1.0);
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<> > rng(eng, normal);
    long N= 1e5;
    long M=opt.T* 50;
    
    sampleCalculator_Euro lsEuro(opt, BS(opt).price());
    SimuLS_CV LSAmerCV(opt, N, M, lsEuro);
    
    std::cout<<LSAmerCV.valuation_calibrated() << std::endl<< LSAmerCV.valuation_stdiv_calibrated()<< std::endl << std::endl<< LSAmerCV.valuation_raw()<<std::endl<< LSAmerCV.valuation_stdiv_raw()<<std::endl;
    
    return 0;

}
 
*/


#endif /* SimuLS_hpp */
