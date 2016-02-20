//
//  SimuNonPathDepEuroBasket.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef SimuNonPathDepEuroBasket_hpp
#define SimuNonPathDepEuroBasket_hpp

#include <stdio.h>
#include "options.hpp"

class SimuNonPathDepEuroBasket {
    // SimuNonPathDepEuroBasket prices a non path dependent basket option by simulations.
    // Since the payoff only relates to asset prices at final stage, we only simulate the final asset prices.
protected:
    long N;
    nonPathDependentBasket_option opt;
    Eigen::MatrixXd asset_price;
    // asset_price stores the asset prices at expiration. It has dimension N*p, where p= opt.count_assets;
    Eigen::VectorXd option_value;
    // option_value has dimension N. It stores the option value at each of the N paths.
    long double mean;
    long double stdiv;
    
public:
    SimuNonPathDepEuroBasket(){};
    SimuNonPathDepEuroBasket(const nonPathDependentBasket_option & o, long paths, unsigned int seed= int(time(0)));
    // With no argument specifying the random number generator to be used, we use boost::mt19937 as random number engine with given seed which has default value transformed from time(0)
    SimuNonPathDepEuroBasket(const nonPathDependentBasket_option & o, long paths, const  Eigen::MatrixXd & RN);
    // RN is a boost matrix containing the CORRELATED std normal random number to be used in the simulation. It must have dimension N*p, p= opt.count_assets;
    
    virtual ~SimuNonPathDepEuroBasket(){};
    
    Eigen::MatrixXd assetPriceDist(){return asset_price;};
    Eigen::VectorXd optionValueDist(){return option_value;};
    
    double valuation(){return mean;};
    double valuation_stdiv(){return stdiv;};
    
};



/*
 
 ******** test *************
 
 class prod_payoff: public ROOT::Math::IBaseFunctionMultiDim{
 protected:
 double K;
 int asset_count;
 bool Call;
 
 double DoEval(const double * x) const{
 double res=1;
 for (int i=0; i<asset_count; i++) res*= x[i];
 
 return fmax(Call?res-K: K-res, 0.0);
 };
 public:
 prod_payoff(){};
 prod_payoff(double strike, int count, bool CallORPut): K(strike), asset_count(count), Call(CallORPut){};
 virtual ~ prod_payoff(){};
 
 unsigned int NDim() const {return asset_count;};
 
 
 
 ROOT::Math::IBaseFunctionMultiDim* Clone() const {
 return new prod_payoff(*this);
 };
 
 };
 
 
 int main(){
 asset a1(1,0.03, 0.01, 0.1);
 asset a2= a1;
 asset a3= a1;
 
 Eigen::MatrixXd Omega(3,3);
 
 Omega(0,0)= 0.01; Omega(0,1)= 0.01*0.6; Omega(0,2)= 0.01*-0.3;
 Omega(1,1)= 0.01; Omega(1,2)=0.01*0.5;
 Omega(2,2)= 0.01;
 
 for (long i=0; i<3; i++)
 for (int j=0; j<i; j++)
 Omega(i,j)= Omega(j,i);
 
 
 prod_payoff payoff(1.2, 3, 1);
 
 std::vector<asset> asset_vec;
 asset_vec.push_back(a1); asset_vec.push_back(a2); asset_vec.push_back(a3);
 
 nonPathDependentBasket_option opt(asset_vec, 1, &payoff, Omega);
 
 SimuNonPathDepEuroBasket simu_EuroBasket(opt, 1e5);
 
 std::cout<< simu_EuroBasket.valuation()<<std::endl<<simu_EuroBasket.valuation_stdiv()<<std::endl;
 
 
 
 return 0;
 };


*/





#endif /* SimuNonPathDepEuroBasket_hpp */
