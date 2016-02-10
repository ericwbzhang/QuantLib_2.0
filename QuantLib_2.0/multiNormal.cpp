//
//  multiNormal.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "multiNormal.hpp"
#include "cholesky.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

multiNormal::multiNormal(boost::numeric::ublas::vector<double> m, boost::numeric::ublas::matrix<double> sigma) {
    mean= m;
    cov= sigma;
    cholesky chol(cov);
    U= chol.factor();
}

boost::numeric::ublas::matrix<double> multiNormal::sample(long n, unsigned int seed) {
    
    long p= mean.size();
    boost::numeric::ublas::matrix<double> X(n, p);
    
    boost::mt19937 eng(seed);
    boost::normal_distribution<double> normal(0.0, 1.0);
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>> rng(eng, normal);
    
    boost::numeric::ublas::vector<double > y(p);
    for (long i= 0; i< n; i++) {
        y.clear();
        for (long j= 0; j< p ; j++) y(j)= rng();
        
        boost::numeric::ublas::row(X, i)= boost::numeric::ublas::prod(boost::numeric::ublas::trans(U), y);
    }
    return X;
}


boost::numeric::ublas::matrix<double> multiNormal::sample(long n, boost::numeric::ublas::matrix<double> RN){
    
    boost::numeric::ublas::matrix<double> X=RN;
    boost::numeric::ublas::vector<double> mr;
    for (long i= 0; i< n; i++) {
       mr= boost::numeric::ublas::row(X, i);
       boost::numeric::ublas::row(X, i)= boost::numeric::ublas::prod(boost::numeric::ublas::trans(U), mr);
    }
    return X;
    
}

