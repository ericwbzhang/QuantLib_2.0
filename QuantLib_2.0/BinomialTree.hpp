//
//  BinomialTree.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 1/28/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef BinomialTree_hpp
#define BinomialTree_hpp

#include <stdio.h>
#include "options_info.hpp"
#include <Eigen/Dense>


class BinomialTree{
protected:
    long N; // number of time steps
    option opt; // the option to be priced
    Eigen::VectorXd value; // the boost vector to capture the option value evolution
    double t, p, q, u, d, disc;
    // t= length of each time step;
    // p= the risk neutral prob of up; q=1-p, the risk neutral prob of down
    // u= the risk neutral multiplier of up; d the risk neutral multiplier of down
    // disc= the risk neutral discounting factor for one time step

public:
    BinomialTree(){};
    BinomialTree( const option &o, long steps);
    
    virtual ~BinomialTree(){};
    
    virtual double valuation();
    virtual double delta();
    virtual double gamma();
    virtual double theta();

};

#endif /* BinomialTree_hpp */
