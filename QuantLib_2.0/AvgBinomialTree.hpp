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
namespace QLib{
    namespace Tree{
        using namespace QLib;
        using namespace QLib::Tree;

        
class AvgBinomialTree{

protected:

    BinomialTree *bt1;
    BinomialTree *bt2;
    
public:
    AvgBinomialTree(){};
    AvgBinomialTree(const option & o, long steps, int tree_type) {
        if (tree_type==0) {
            // BinomialTree
            bt1= new BinomialTree(o, steps);
            bt2= new BinomialTree(o, steps -1);
        }else {
            // BinomialBSTree
            bt1= new BinomialBSTree(o, steps);
            bt2= new BinomialBSTree(o, steps-1);
        }
        

    };
    // tree_type=0 if we want average BinomialTree;
    // tree_type=1 if we want average BinomialBSTree;
    
    virtual ~AvgBinomialTree() {delete bt1; delete bt2; };
    
    double valuation(){return .5*(bt1->valuation()+ bt2-> valuation());};
    double delta(){return .5*(bt1->delta()+ bt2->delta());};
    double gamma(){return .5*(bt1->gamma()+ bt2-> gamma());};
    double theta(){return .5*(bt1->theta()+ bt2->theta());}
;
    
};
    }
}

#endif /* AvgBinomialTree_hpp */
