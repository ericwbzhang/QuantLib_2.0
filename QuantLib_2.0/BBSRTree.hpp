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
namespace QLib{
    namespace Tree{
        using namespace QLib;
        using namespace QLib::Tree;

        
class BBSRTree {
protected:
    BinomialBSTree *bt1;
    BinomialBSTree *bt2;
    // BBSRTree method involves two BinomialBSTree, one has steps N, another has steps N/2. Orginial, BinomialBSTree converges as O(1/N), while BBSRTree converges as O(1/N^2)
    
public:
    BBSRTree(){};
    BBSRTree(const option & o, long steps){
        bt1= new BinomialBSTree(o, steps);
        bt2= new BinomialBSTree(o, long(steps/2));
    };
    virtual ~BBSRTree(){delete bt1; delete bt2;} ;
    
    double valuation(){    return 2*bt1->valuation()- bt2->valuation();};
    double delta(){
        return  2*bt1->delta()- bt2->delta();
    };
    double gamma(){
        return 2*bt1->gamma()- bt2->gamma();
    };
    double theta(){return 2*bt1->theta()- bt2->theta();};
    
};
    }
}

#endif /* BBSRTree_hpp */
