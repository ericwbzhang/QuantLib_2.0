//
//  test.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include <iostream>
#include "miscellaneous.h"
#include "option_TreePricers.h"
#include "BS.hpp"
#include "option_Simulation.h"
#include <boost/random.hpp>
#include "algebra.h"

class prod_payoff: public realValueFunctor{
protected:
    double K;
    
public:
    prod_payoff(){};
    prod_payoff(double strike): K(strike){};
    
    virtual double operator()(const std::vector<double>&  args){
        double res=1;
        for (std::vector<const double>::iterator it= args.begin(); it!= args.end(); it++){
            res*= *it;
        }
        
        return fmax(0.0, res-K);
    };
    
    virtual double operator()(const Eigen::VectorXd & args){
        long n=args.size();
        std::vector<double> vec(n);
        for (long i=0; i<n; i++) vec[i]= args(i);
        
        return this-> operator()(vec);
    };
    
    virtual realValueFunctor* clone(){
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
    
    
    prod_payoff payoff(1.2);
    
    std::vector<asset> asset_vec;
    asset_vec.push_back(a1); asset_vec.push_back(a2); asset_vec.push_back(a3);
    
    nonPathDependentBasket_option opt(asset_vec, 1, &payoff, Omega);
    
    SimuNonPathDepEuroBasket simu_EuroBasket(opt, 1e5);
    
    std::cout<< simu_EuroBasket.valuation()<<std::endl<<simu_EuroBasket.valuation_stdiv()<<std::endl;

    
    std::srand(int(time(0)));
    Eigen::MatrixXd A(4,4);
    A<< 10.,-1.,2., 0.,
        -1., 11., -1., 3.,
        2., -1., 10., -1.,
        0.0, 3., -1., 8.;
    Eigen::MatrixXd B(4,1);
    B<< 6., 25., -11., 15.; 
    linearSolver_SOR sol(A, B, 1e-8);
    

  
    std::cout<<sol.iterations()<<std::endl;
    std::cout<< sol.solve()<< std::endl<<std::endl;
    std::cout<< A*sol.solve()-B<<std::endl;
//    std::cout<< qr.matrixQ()<<std::endl<<std::endl;
//    std::cout<< qr.matrixR()<< std::endl<<std::endl;
  //  std::cout<< (svd.matrixU()*svd.matrixS()* (svd.matrixV().transpose())- X).lpNorm<Eigen::Infinity>()<<std::endl;
    
    
    //std::cout<< lu_f.determinant()<<std::endl;
    //  if we are taking determinant wrt non square matrix, above will show runtime error. 
    return 0;
    
    
    
};
