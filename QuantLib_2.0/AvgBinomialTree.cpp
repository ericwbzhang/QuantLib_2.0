//
//  AvgBinomialTree.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "AvgBinomialTree.hpp"

AvgBinomialTree::AvgBinomialTree(option o, long steps, int tree_type){
    if (tree_type==0) {
        // BinomialTree
        bt1= new BinomialTree(o, steps);
        bt2= new BinomialTree(o, steps -1);
    }else {
        // BinomialBSTree
        bt1= new BinomialBSTree(o, steps);
        bt2= new BinomialBSTree(o, steps-1);
    }
    
}

double AvgBinomialTree::valuation() {return .5*(bt1->valuation()+ bt2-> valuation());}


double AvgBinomialTree::delta() {return .5*(bt1->delta()+ bt2->delta());}

double AvgBinomialTree::gamma() {return .5*(bt1->gamma()+ bt2-> gamma());}

double AvgBinomialTree::theta() {return .5*(bt1->theta()+ bt2->theta());}

