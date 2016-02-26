//
//  FD_vanilla.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/25/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef FD_vanilla_hpp
#define FD_vanilla_hpp

#include <stdio.h>
#include "HeatPDE.hpp"
#include "options_info.hpp"

namespace QLib {
    namespace FD{

class FD_vanilla {
// FD_vanilla price vanilla options (European and American) using finite difference method
    

protected:
    // ****** data member
    long N, M;
    double a, b; // N, M ,a ,b are params to use in computation
    double SMIN, SMAX , tMIN; // SMIN, SMAX and tMIN defines the valuation horizon on (S, t) space
    double xMIN, xMAX, tauMAX;  // xMAX, xMIN, tauMAX are on the (x, tau) space,  correpond to SMAX, SMIN and tMIN.
    option opt;
    // tau0Margin, xLeftMargin, xRightMargin are the boundary conditions on (x, tau) space
    std::shared_ptr<ROOT::Math::IBaseFunctionOneDim> tau0Margin;
    std::shared_ptr<ROOT::Math::IBaseFunctionOneDim> xLeftMargin;
    std::shared_ptr<ROOT::Math::IBaseFunctionOneDim> xRightMargin;
    std::shared_ptr<ROOT::Math::IBaseFunctionMultiDim> executionBenefit;
    std::shared_ptr<HeatPDE> pdeKernal;
    
    //***** method
    virtual void Euro    (int algorithm, double weight, double tol, long iterMAX){
        pdeKernal.reset(new HeatPDE(1, xMIN, xMAX, tauMAX,
                                    N, M,
                                    * tau0Margin,
                                    * xLeftMargin,
                                    * xRightMargin,
                                    algorithm, weight, tol, iterMAX));
        

    };
    virtual void American(int algorithm, double weight, double tol, long iterMAX){
        pdeKernal.reset(new HeatPDE_gridUpdate(1, xMIN, xMAX, tauMAX,
                                               N, M,
                                               *tau0Margin,
                                               *xLeftMargin,
                                               *xRightMargin,
                                               * executionBenefit,
                                               algorithm, weight, tol, iterMAX));
        

    };
    virtual void transf_St2xtau(double S, double t, double *xtau) const {
        xtau[0]= log(S/opt.K);
        xtau[1]= opt.sigma* opt.sigma*.5 * (opt.T- t);
    };
    virtual void transf_xtau2St(double x, double tau, double *St) const {
        St[0]= opt.K* exp(x);
        St[1]= opt.T- 2*tau/(opt.sigma*opt.sigma);
    };
    virtual double transf_u2value(double u, double x, double tau)const {
        return exp(-a*x-b*tau)*u;
    };
    

    // subclass
    
    class exeBenefit: public ROOT::Math::IBaseFunctionMultiDim {
    protected:
        bool Call;
        double K, A, B;
        
        double DoEval (const double * y) const {
            return Call?
            fmax(K* exp(y[0])- K, 0.0 )* exp(A*y[0]+ B* y[1]):
            fmax(K- K* exp(y[0]), 0.0 )* exp(A*y[0]+ B* y[1]);
        };
    public:
        exeBenefit (double k, double a, double b, bool call): K(k), A(a), B(b), Call(call){};
        unsigned int NDim() const {return 2;};
        ROOT::Math::IBaseFunctionMultiDim * Clone () const {return new exeBenefit(*this);};
    };
    
    class f: public ROOT::Math::IBaseFunctionOneDim{
    
    protected:
        double A;
        double K;
        bool Call;
        
        double DoEval(double x) const {
            return exp(A*x)*
            (Call?  fmax(exp(x)*K -K, 0.0 ):
                    fmax(K- exp(x)*K, 0.0));
        };
    public:
        f( double a, double k, bool C) : A(a),K(k), Call(C){};
        
        ROOT::Math::IBaseFunctionOneDim * Clone()const{
            return new f(*this);
        };
    };
    
    class g_left: public ROOT::Math::IBaseFunctionOneDim{
    protected:
        bool Call;
        double C1, C2, C3, C4;
      
        double DoEval(double x) const {
            return Call? 0:
            C1* exp(x*C2)- C3* exp(x*C4);
        };
    public:
        g_left(double * C, bool call): Call(call), C1(C[0]), C2(C[1]), C3(C[2]), C4(C[3]){};
        
        ROOT::Math::IBaseFunctionOneDim * Clone() const{return new g_left(*this);};
    };
    
    
    class g_right :public ROOT::Math::IBaseFunctionOneDim{
    protected:
        bool Call;
        double C1, C2, C3,  C4;
        
        double DoEval(double x) const {
            return Call?
            C1*exp(x*C2)- C3* exp(x*C4):
            0.0 ;
        };
    public:
        g_right(double * C,  bool call): Call(call), C1(C[0]), C2(C[1]), C3(C[2]), C4(C[3]){};
    
        ROOT::Math::IBaseFunctionOneDim* Clone() const{return  new g_right(*this);};
    };
    
public:
    
    //constructor
    FD_vanilla(){};
    // use default copy constructor
    FD_vanilla(const option & o, double priceMin, double priceMax, double tMin,
               long n, long m, int algorithm,
               double weight= 1.1, double tol=1.e-7, long iterMAX= 50){
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

    };
    // destructor
    virtual ~ FD_vanilla(){};
    
    // methods
    
    double valuation(double S, double t) const{
        double xtau[2];
        this->transf_St2xtau(S,t, xtau);
        
        return this->transf_u2value(pdeKernal->evaluation(xtau[0], xtau[1]), xtau[0], xtau[1]);

    };
    double delta(double S, double t) const{
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

    };
    double gamma(double S, double t) const{
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
        

    };
    double theta(double S, double t) const{
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
        

    };
    
    int computationMethod() const   {return pdeKernal->computationMethod();};
    double computationAlpha() const {return pdeKernal-> computationAlpha();};
    bool computationSuccess() const {return pdeKernal-> computationSuccess();};
    
    
    long logPriceDiscret() const {return N;};
    long expirationDiscret() const {return M;};
    option optionInfo() const {return opt;};
};






/*
 
 // ************* test ********

int main(){
    
    option opt(44,40, 1, 0.06, 0, 0.4, 1, 1); // S=K=1, T=1, r=q=0.01, sigma= 0.2, Call=1, Euro=1
    
    opt.Call =0;
    
    opt.Euro=0;
    opt.T=2;
    
    option_BS bs(opt);
    std:: cout<< bs.price()<< std::endl;
    

    
    double pMax= opt.S* exp(5* opt.sigma* sqrt(opt.T));
    double pMin= opt.S* exp(-5* opt.sigma*  sqrt(opt.T));
    long n= 1e2;
    long m= 1e4;
    FD_vanilla fd(opt, pMin, pMax, 0, n, m, 2);
    
    std::cout<< fd.computationAlpha()<< std::endl<<fd.computationSuccess()<<std::endl<<fd.valuation(opt.S, 0)<<std::endl;
    
    
    return 0;
    
};

*/
    }
}

#endif /* FD_vanilla_hpp */
