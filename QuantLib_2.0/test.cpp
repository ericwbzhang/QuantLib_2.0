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
#include "option_FD.h"
#include <boost/random.hpp>
#include "algebra.h"
#include "lm.hpp"
#include <random>
#include <Math/IFunction.h> // ROOT Math library

using namespace QLib;
class func_f: public ROOT::Math::IBaseFunctionOneDim {
protected:
    double DoEval(double x) const {
        return 1-x*x;
    };
    
public:
    ROOT::Math::IBaseFunctionOneDim * Clone() const {
        return new func_f(*this);
    };
};

class func_g :public ROOT::Math::IBaseFunctionOneDim {
protected:
    double DoEval(double x) const{
        return 0;
    };
    
public:
    ROOT::Math::IBaseFunctionOneDim * Clone() const{return  new func_g(*this);};
};

int main(){
        
    option opt(44,40, 1, 0.06, 0, 0.4, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
    
    opt.Call =0;
    
//    
//    func_f f;
//    func_g g;
//    HeatPDE heatpde(1, -1, 1, 2, 1e2, 1e4, f, g, g, 1);
//    
//    std::shared_ptr<ROOT::Math::IBaseFunctionMultiDim> heatpde_formula;
//    Eigen::MatrixXd grid= heatpde.computationGrid();
//    
//    
//    std::cout<< heatpde.computationAlpha()<< std::endl<< heatpde.computationSuccess()<< std::endl << std::endl;
//    
//    std::cout<< heatpde.evaluation(0, .1)<<std::endl;
//    

    
    opt.Euro=0;
    opt.T=2;
    
    option_BS bs(opt);
    std:: cout<< bs.price()<< std::endl;

    long N= 1e5;
    long M=opt.T* 50;

    double pMax= opt.S* exp(5* opt.sigma* sqrt(opt.T));
    double pMin= opt.S* exp(-5* opt.sigma*  sqrt(opt.T));
    long n= 1e2;
    long m= 1e4;
    QLib::FD::FD_vanilla fd(opt, pMin, pMax, 0, n, m, 2);
    
    std::cout<< fd.computationAlpha()<< std::endl<<fd.computationSuccess()<<std::endl<<fd.valuation(opt.S, 0)<<std::endl;
    
    Simu::sampleCalculator_Euro lsEuro(opt, option_BS(opt).price());
    Simu::SimuLS_CV LSAmerCV(opt, N, M, lsEuro);
    
    std::cout<<LSAmerCV.valuation_calibrated() << std::endl<< LSAmerCV.valuation_stdiv_calibrated()<< std::endl << std::endl<< LSAmerCV.valuation_raw()<<std::endl<< LSAmerCV.valuation_stdiv_raw()<<std::endl;
    
    Tree::BBSRTree (opt, 1e4);
    


    return 0;
    
};
