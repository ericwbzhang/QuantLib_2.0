//
//  QR_HH.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/12/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef QR_HH_hpp
#define QR_HH_hpp

#include <stdio.h>
#include <Eigen/Dense>
namespace QLib{
    using namespace QLib;

class QR_HH{
// QR_HH wraps Eigen::HouseholderQR. It takes input matrix A which is a general real matrix.
// For general real matrix A, it has a decomp A= QR. Q is an orthogonal matrix and R is an upper triangular matrix. If A is invertible and we require the diag of A to be positive, the decomp is unique.
protected:
    Eigen::HouseholderQR<Eigen::MatrixXd> qr;
    
public:
    
    QR_HH(){};
    QR_HH(const Eigen::MatrixXd &A) {qr= A.householderQr();};
    virtual ~QR_HH(){};
    
    Eigen::MatrixXd matrixQ(){
        return qr.householderQ();
    };
    
    Eigen::MatrixXd matrixR(){
        return Eigen::MatrixXd (qr.matrixQR(). triangularView<Eigen::Upper>());
    };
    
    
    
};
}
#endif /* QR_HH_hpp */
