//
//  AvgBinomialTree.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef AvgBinomialTree_hpp
#define AvgBinomialTree_hpp

#include <stdio.h>
#include "BinomialBSTree.hpp"

class AvgBinomialTree{

protected:

    BinomialTree *bt1;
    BinomialTree *bt2;
    
public:
    AvgBinomialTree(){};
    AvgBinomialTree(option o, long steps, int tree_type);
    // tree_type=0 if we want average BinomialTree;
    // tree_type=1 if we want average BinomialBSTree;
    
    virtual ~AvgBinomialTree() {delete bt1; delete bt2; };
    
    double valuation();
    double delta();
    double gamma();
    double theta();
    
};

#endif /* AvgBinomialTree_hpp */
