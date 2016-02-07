//
//  BBSRTree.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef BBSRTree_hpp
#define BBSRTree_hpp

#include <stdio.h>

#include "BinomialBSTree.hpp"

class BBSRTree {
protected:
    BinomialBSTree *bt1;
    BinomialBSTree *bt2;
    // BBSRTree method involves two BinomialBSTree, one has steps N, another has steps N/2. Orginial, BinomialBSTree converges as O(1/N), while BBSRTree converges as O(1/N^2)
    
public:
    BBSRTree(){};
    BBSRTree(option o, long steps);
    virtual ~BBSRTree(){delete bt1; delete bt2;} ;
    
    double valuation();
    double delta();
    double gamma();
    double theta();
    
};

#endif /* BBSRTree_hpp */
