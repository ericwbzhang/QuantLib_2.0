//
//  QR_HHColPiv.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef QR_HHColPiv_hpp
#define QR_HHColPiv_hpp

#include <stdio.h>
#include <Eigen/Dense>

namespace QLib{
    using namespace QLib;

class QR_HHColPiv{
// QR_HHColPiv wrapps Eigen::ColPivHouseHolderQR. It takes input matrix A, a general real (or complex) matrix, and decompose it as AP= QR. P is a permutation matrix (we do col pivoting), and Q is an orthognal (unitary in complex) matrix and R is an upper triangular matrix.
// QR decomp applies to general matrix. The col pivoting QR generates more stable result and higher accuracy.
    
protected:
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_c;
    
public:
    QR_HHColPiv(){};
    QR_HHColPiv(const Eigen::MatrixXd &A) {
        qr_c= A.colPivHouseholderQr();
    };
    
    virtual ~QR_HHColPiv(){};
    
    Eigen::PermutationMatrix<Eigen::Dynamic> matrixP() {return qr_c.colsPermutation();};
    Eigen::MatrixXd matrixQ(){
        Eigen::MatrixXd Q=qr_c.householderQ();
        return Q;
    };
    
    Eigen::MatrixXd matrixR(){
        Eigen::MatrixXd R(qr_c.matrixQR().triangularView<Eigen::Upper>());
        return R;
    };
    
    long rank() { return qr_c.rank();};
    
    bool isInvertible () {return qr_c.isInvertible();};
    
    Eigen::MatrixXd inverse(){
        // If A is not invertible, the return will have undefined entries or error. Check whether invertible before call the inverse method. 
        return qr_c.inverse();};
};


}
#endif /* QR_HHColPiv_hpp */
