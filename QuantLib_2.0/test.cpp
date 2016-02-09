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

int main(int argc, const char * argv[]) {
    

    option opt(1,1, 1, 0.01, 0.01, 0.2, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
    barrier_option b_opt(opt, 1.1, 1);

    return 0;
}
 