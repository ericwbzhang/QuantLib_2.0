//
//  FDEEuro.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/21/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef FDEEuro_hpp
#define FDEEuro_hpp

#include <stdio.h>
#include "HeatPDE.hpp"
#include "options_info.hpp"

class FDEEuro {
protected:
    option opt;
    
    HeatPDE PDEsolver;
    
public:
    FDEEuro
};

#endif /* FDEEuro_hpp */
