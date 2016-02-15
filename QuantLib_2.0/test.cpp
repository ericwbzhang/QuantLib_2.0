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
#include "lm.hpp"
#include <dlib/array.h>
#include <random> 


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
    
//    Eigen::MatrixXd tmp( LSAmer.optionCashFlow());
//    for (long i= 0; i< tmp.rows(); i++){
//        int count=0;
//        for (long j=0; j<tmp.cols(); j++){
//            if (tmp(i,j)) count++;
//        }
//        if (count>1) flag= false;
//    }
//    
//    std::cout<<flag<<std::endl;
    
//    std::srand(int(time(0)));
//    Eigen::MatrixXd A;
//    A.setRandom(1000, 1)*1e4;
//    // Eigen::MatrixXd C= A+ A.transpose();
//    Eigen::MatrixXd B;
//    B.setRandom(1000,1);B= 2*A+B;
//    Eigen::MatrixXd C(1000,2); C<< A, A;
//    C.col(0).setOnes();
//
//    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
//    std::cout<< qr.solve(B)<<std::endl<<std::endl;
//    
//    lm_Ridge lm(C,B);
//    std::cout<< lm<<std::endl;
      // A = perm * A; // permute rows
    //std::cout<<A;
//
//    std::cout<<es.isEigenValueReal()<<std::endl;
//    std::cout<< es.eigenvalues_real()<< std::endl<<std::endl;
//    std::cout<< es.eigenvectors_real()* es.eigenvalues_real().asDiagonal()* es.eigenvectors_real().inverse()- A<<std::endl;
//    std::cout<< qr.matrixQ()<<std::endl<<std::endl;
//    std::cout<< qr.matrixR()<< std::endl<<std::endl;
  //  std::cout<< (svd.matrixU()*svd.matrixS()* (svd.matrixV().transpose())- X).lpNorm<Eigen::Infinity>()<<std::endl;
    
    
    //std::cout<< lu_f.determinant()<<std::endl;
    //  if we are taking determinant wrt non square matrix, above will show runtime error. 
    return 0;
    
};
