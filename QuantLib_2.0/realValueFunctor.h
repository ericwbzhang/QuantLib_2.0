//
//  realValueFunctor.h
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef realValueFunctor_h
#define realValueFunctor_h

#include <boost/numeric/ublas/vector.hpp> 
#include <stdio.h>

class realValueFunctor{
    // realValueFunctor is a virtual class. It sets the interface for general function from R^p to R.

public:
    realValueFunctor(){};
    virtual ~realValueFunctor(){};
    
    virtual double operator() (std::vector<double> args) =0;
    virtual double operator() (boost::numeric::ublas::vector<double> args )=0;
    
    virtual realValueFunctor* clone()  =0;
};




#endif /* realValueFunctor_h */
