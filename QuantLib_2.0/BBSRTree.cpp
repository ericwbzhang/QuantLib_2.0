//
//  BBSRTree.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/7/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "BBSRTree.hpp"

BBSRTree::BBSRTree(option o, long steps) {
    bt1= new BinomialBSTree(o, steps);
    bt2= new BinomialBSTree(o, long(steps/2));
    
}

double BBSRTree::valuation() {
    return 2*bt1->valuation()- bt2->valuation();
}

double BBSRTree::delta() {
    return  2*bt1->delta()- bt2->delta();
}

double BBSRTree:: gamma() {
    return 2*bt1->gamma()- bt2->gamma();
}

double BBSRTree::theta(){return 2*bt1->theta()- bt2->theta();}