//
//  test.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include <iostream>
#include "option_TreePricers.h"
#include "option_Simulation.h"
#include <boost/random.hpp>
#include "algebra.h"
#include "lm.hpp"
#include <random>
#include <Math/IFunction.h> // ROOT Math library

class func: public ROOT::Math::IBaseFunctionOneDim{

    double DoEval(double x) const {return 1;};
    
    ROOT::Math::IBaseFunctionOneDim * Clone() const {return new func(); };
};

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
        
//    option opt(36,40, 1, 0.06, 0, 0.4, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
//    
//    opt.Call =0;
//    option_BS bs(opt);
//    std:: cout<< bs.price()<< std::endl;
//    
//    opt.Euro=0;
//    
//    boost::mt19937 eng(time(0));
//    boost::normal_distribution<> normal(0,1.0);
//    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<> > rng(eng, normal);
//    long N= 1e5;
//    long M=opt.T* 50;
//
//    sampleCalculator_Euro lsEuro(opt, option_BS(opt).price());
//    SimuLS_CV LSAmerCV(opt, N, M, lsEuro);
//    
//    std::cout<<LSAmerCV.valuation_calibrated() << std::endl<< LSAmerCV.valuation_stdiv_calibrated()<< std::endl << std::endl<< LSAmerCV.valuation_raw()<<std::endl<< LSAmerCV.valuation_stdiv_raw()<<std::endl;
//    

    
        std::srand(int(time(0)));
        Eigen::MatrixXd A;
        A.setRandom(1000, 1)*1e4;
        // Eigen::MatrixXd C= A+ A.transpose();
        Eigen::MatrixXd B;
        B.setRandom(1000,1);B= 2*A+B;
        Eigen::MatrixXd C(1000,2); C<< A, A;
        C.col(0).setOnes();
    
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        std::cout<< qr.solve(B)<<std::endl<<std::endl;
    
        lm_Ridge lm(C,B);
        std::cout<< lm<<std::endl;
    
    
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
