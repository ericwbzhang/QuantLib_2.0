//
//  cholesky.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "cholesky.hpp"

cholesky::cholesky(boost::numeric::ublas::matrix<double> A) {
    // A is a spd matrix with size n*n
    
    long n= A.size1();
    U.resize(n, n);
    U.clear();
    
    if (n== A.size2()){
        flag= false;
        
        for (long i=0; i<n-1; i++) {
            U(i,i)= sqrt(A(i,i));
            flag= flag|| isnan(U(i,i));
            
            for (long k= i+1; k< n; k++) {
                U(i,k)= A(i,k)/ U(i,i);
                flag= flag || isnan(U(i,k));
            }
            
            for (long j=i+1; j<n; j++) {
                for (long k=j; k<n; k++){
                    A(j,k)-=U(i,j)*U(i,k);
                }
            }
        }
        
        U(n-1, n-1)= sqrt(A(n-1, n-1));
        flag=flag|| isnan(U(n-1, n-1));
        
    }else {
        // A is not a square matrix;
        flag= true;
    }
}