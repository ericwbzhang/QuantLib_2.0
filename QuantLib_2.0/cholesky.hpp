//
//  cholesky.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/9/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef cholesky_hpp
#define cholesky_hpp

#include <stdio.h>
#include <boost/numeric/ublas/matrix.hpp>

class cholesky{
// cholesky plays cholesky decomposition to the input matrix A, which mush be symmetric positive definite, and return the factor matrix U, which is upper triangle. A= U'U;

protected:
    boost::numeric::ublas::matrix<double> U;
    bool flag; // flag= true if there is error appearing in the decomposition, which means the input is not spd.
    
public:
    cholesky(){};
    cholesky(boost::numeric::ublas::matrix<double> A);
    
    virtual ~cholesky(){};
    
    boost::numeric::ublas::matrix<double> factor() {return U; };
    bool spd(){return !flag;};

};

/*
 ********** test***********

int main(){

    boost::numeric::ublas:: matrix<double> A(4,4);
    A(0,0)= 9; A(0,1)= -3; A(0,2)=6; A(0,3)= -3;
    A(1,1)=5; A(1,2)= -4; A(1,3)= 7;
    A(2,2)=21; A(2,3)= 3;
    A(3,3)= 15;
    
    for (long i= 0 ;i <4; i++) {
        for (int j= 0; j<i; j++)
            A(i,j)= A(j,i);
    }
    
    cholesky chol(A);
    
    boost::numeric::ublas::matrix<double> U= chol.factor();
    boost::numeric::ublas::matrix<double> X= boost::numeric::ublas::trans(U);
    X= boost::numeric::ublas::prod(X, U);
    
    std::cout<< boost::numeric::ublas::norm_1(A- X)<<std::endl;
    
    return 0;
}
*/


#endif /* cholesky_hpp */
