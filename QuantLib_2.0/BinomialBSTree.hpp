//
//  BinomialBSTree.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef BinomialBSTree_hpp
#define BinomialBSTree_hpp

#include <stdio.h>
#include "BinomialTree.hpp" 
#include "BS.hpp"

class BinomialBSTree: public BinomialTree {
// The content of Binomial Tree
/*
protected:
    long N; // number of time steps
    option opt; // the option to be priced
    boost::numeric::ublas::vector<double> value; // the boost vector to capture the option value evolution
    double t, p, q, u, d, disc;
    // t= length of each time step;
    // p= the risk neutral prob of up; q=1-p, the risk neutral prob of down
    // u= the risk neutral multiplier of up; d the risk neutral multiplier of down
    // disc= the risk neutral discounting factor for one time step
*/
    
public:
    BinomialBSTree(){};
    BinomialBSTree(option o, long steps ) ;
    
    virtual ~BinomialBSTree(){};
    
    
};


#endif /* BinomialBSTree_hpp */
