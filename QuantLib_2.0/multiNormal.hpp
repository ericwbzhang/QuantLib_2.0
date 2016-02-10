//
//  multiNormal.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef multiNormal_hpp
#define multiNormal_hpp

#include <stdio.h>
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>

class multiNormal{
    // multiNormal generates multivariate normal distributed RN given the mean vector and covariance matrix. It utlizes the one dimension RN generator in boost (mt19937 random engine).
protected:
    boost::numeric::ublas::vector<double> mean;
    boost::numeric::ublas::matrix<double> cov;
    boost::numeric::ublas::matrix<double> U;
    // cov= U'U;
    
public:
    multiNormal(){};
    multiNormal( boost::numeric::ublas::vector<double> m, boost::numeric::ublas::matrix<double> sigma);
    
    virtual ~multiNormal(){};
    
    boost::numeric::ublas::matrix<double> sample(long n, unsigned int seed= int(time(0)));
    // it returns a sample of size n. Each record has p elements, p= cov.size1();
    
    boost::numeric::ublas::matrix<double> sample(long n, boost::numeric::ublas::matrix<double> RN);
    // RN is a matrix of n*p. It contains iid std normal distributed RNs. 
    
};

#endif /* multiNormal_hpp */
