//
//  FD_vanilla.cpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/25/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#include "FD_vanilla.hpp"

FD_vanilla::FD_vanilla(const option & o, double priceMin, double priceMax, double tMin, long n, long m, int algorithm, double weight, double tol, long iterMAX){

    opt= o;
    N=n;
    M=m;
    SMAX= priceMax;
    SMIN= priceMin;
    tMIN= tMin;
    a= (opt.r- opt.q)/ (opt.sigma* opt.sigma) -.5;
    b= 2* opt.q/opt.sigma +pow((opt.r-opt.q)/ (opt.sigma* opt.sigma)+.5, 2.0);
    
    double CL[4];
    double CR[4];
    xMIN= log(SMIN/ opt.K);
    xMAX= log(SMAX/ opt.K);
    tauMAX= (opt.T- tMIN)/2.0 * opt.sigma* opt.sigma;

    
    CL[0]= opt.K* exp(a* xMIN);
    CL[1]= b- 2* opt.r/ (opt.sigma* opt.sigma);
    CL[2]= CL[0]* exp(xMIN);
    CL[3]= b- 2* opt.q/ (opt.sigma* opt.sigma);
    
    CR[1]= CL[3];
    CR[2]= opt.K* exp(a* xMAX);
    CR[3]= CL[1];
    CR[0]= CR[2]* exp(xMAX);
    
    tau0Margin.     reset(new f(a,opt.K, opt.Call));
    xLeftMargin.    reset(new g_left(CL, opt.Call));
    xRightMargin.   reset(new g_right(CR, opt.Call));
    executionBenefit.reset(new exeBenefit(opt.K,a, b, opt.Call));

    if (opt.Euro) {
        this-> Euro(algorithm, weight, tol, iterMAX);
    }else {
        this->American(algorithm, weight, tol, iterMAX);
    }
    
}

void FD_vanilla::Euro(int algorithm, double weight, double tol, long iterMAX){
    pdeKernal.reset(new HeatPDE(1, xMIN, xMAX, tauMAX,
                            N, M,
                            * tau0Margin,
                            * xLeftMargin,
                            * xRightMargin,
                                algorithm, weight, tol, iterMAX));

}

void FD_vanilla::American(int algorithm, double weight, double tol, long iterMAX){
    pdeKernal.reset(new HeatPDE_gridUpdate(1, xMIN, xMAX, tauMAX,
                                           N, M,
                                           *tau0Margin,
                                           *xLeftMargin,
                                           *xRightMargin,
                                           * executionBenefit,
                                           algorithm, weight, tol, iterMAX));

}


double FD_vanilla::valuation(double S, double t) const {
    double xtau[2];
    this->transf_St2xtau(S,t, xtau);
    
    return this->transf_u2value(pdeKernal->evaluation(xtau[0], xtau[1]), xtau[0], xtau[1]);
}

double FD_vanilla:: delta(double S, double t) const {
    if (S>=SMAX && opt.Call) return 1;
    else if (S>=SMAX && !opt.Call) return 0;
    else if (S<=SMIN && opt.Call) return 0;
    else if (S<=SMIN && !opt.Call) return  -1;
    
    double xtau[2];
    this->transf_St2xtau(S, t, xtau);
    
    double dx= pdeKernal->dx();
    double xtau2[2];
    xtau2[0]=xtau[0]+dx;
    if (xtau2[0]> xMAX ) return opt.Call? 1 :0;
    xtau2[1]=xtau[1];
    
    double v1= this->transf_u2value(pdeKernal->evaluation(xtau[0], xtau[1]), xtau[0], xtau[1]);
    double v2= this->transf_u2value(pdeKernal->evaluation(xtau2[0], xtau2[1]), xtau2[0], xtau2[1]);
    
    double S2[2];
    this->transf_xtau2St(xtau2[0], xtau2[1], S2);
    return (v2-v1)/ (S2[0]- S);

}


double FD_vanilla:: gamma(double S, double t) const {

    double xtau1[2];
    double xtau2[2];
    double xtau3[2];
    
    this->transf_St2xtau(S, t, xtau2);
    double dx= pdeKernal->dx();
    xtau1[0]= xtau2[0]-dx;
    xtau1[1]= xtau2[1];
    xtau3[0]= xtau2[0]+dx;
    xtau3[1]= xtau2[1];
    
    double St1[2];
    double St3[2];
    this->transf_xtau2St(xtau1[0], xtau1[1], St1);
    this->transf_xtau2St(xtau3[0], xtau3[1], St3);
    double delta1= this->delta(St1[0], St1[1]);
    double delta3= this->delta(St3[0], St3[1]);
    
    return (delta3-delta1)/(St3[0]- St1[0]);
    
}

double FD_vanilla:: theta(double S, double t) const {
    // t > tMIN
    double xtau[2];
    this->transf_St2xtau(S, t, xtau);
    
    double dt= pdeKernal->dt();
    double xtau2[2];
    xtau2[1]=xtau[1]+dt;
    xtau2[0]=xtau[0];
    
    double v1= this->transf_u2value(pdeKernal->evaluation(xtau[0], xtau[1]), xtau[0], xtau[1]);
    double v2= this->transf_u2value(pdeKernal->evaluation(xtau2[0], xtau2[1]), xtau2[0], xtau2[1]);
    
    double S2[2];
    this->transf_xtau2St(xtau2[0], xtau2[1], S2);
    return (v2-v1)/ (S2[1]- t);
    
}




