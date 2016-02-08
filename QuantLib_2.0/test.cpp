//
//  test.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include <iostream>
#include "option_TreePricers.h"
#include "BS.hpp"
#include "option_Simulation.h"
#include <boost/random.hpp>

int main(int argc, const char * argv[]) {
    

    option opt(1,1, 1, 0.01, 0.01, 0.2, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
    
    BS bs_formula(opt);
    
    BinomialTree bio_tree(opt, 1e4);
    std::cout<< bs_formula.price()<< std::endl<<"******"<<std::endl<<bio_tree.valuation()<< std::endl;
    
    SimuEuro simu(opt, 1e7);
    
    std::cout<< simu.valuation()<<std::endl<<simu.stdiv()<<std::endl;
    
    return 0;
}
