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


class realValueFunctor{
    // realValueFunctor is a virtual class. It sets the interface for general function from R^p to R.
protected:
    long p; //the dimension of space the function is working on.
public:
    realValueFunctor(){};
    realValueFunctor(long dim): p(dim){};
    virtual ~realValueFunctor(){};
    
    virtual double operator() (std::vector<double> args) =0;
    virtual double operator() (boost::numeric::ublas::vector<double> args )=0;
    virtual long support_dim() {return p; };
};




#endif /* realValueFunctor_h */
