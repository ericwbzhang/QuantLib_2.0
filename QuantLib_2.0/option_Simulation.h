//
//  option_Simulation.h
//  QuantLib_2.0
//
//  Created by EricZ on 2/8/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef option_Simulation_h
#define option_Simulation_h

#include "SimuEuro.hpp"
#include "SimuEuroBarrier.hpp"


/*
 ********* test ********
 
 int main(){
 option opt(1,1, 1, 0.01, 0.01, 0.2, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
 barrier_option b_opt(opt, 1.1, 1);
 
 BS bs_formula(opt);
 std::cout<<bs_formula.price()<< std::endl<< "***********" <<std::endl;
 
 SimuEuro simu(opt, 1e5);
 
 std::cout<< simu.valuation()<<std::endl<<simu.valuation_stdiv()<<std::endl<<"*************"<<std::endl;
 
 
 SimuEuroBarrier simuBarrier(b_opt, 1e5, 2e2);
 std::cout<< simuBarrier.valuation()<<std::endl<< simuBarrier.valuation_stdiv()<<std::endl<<simuBarrier.effectRatio()<<std::endl;
 barrier_option b_opt2= b_opt;
 b_opt2.knock_out=0;
 SimuEuroBarrier simuBarrier2(b_opt2, 1e5, 2e2);
 
 std::cout<< simuBarrier.valuation()+ simuBarrier2.valuation()<<std::endl<< simuBarrier2.valuation_stdiv()+ simuBarrier.valuation_stdiv()<<std::endl;
 return 0;
 
 }
 

 */


#endif /* option_Simulation_h */
