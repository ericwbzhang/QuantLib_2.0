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
    virtual void Euro    (int algorithm, double weight, double tol, long iterMAX);
    virtual void American(int algorithm, double weight, double tol, long iterMAX);
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
    FD_vanilla(const option & o, double priceMin, double priceMax, double tMin, long n, long m, int algorithm, double weight= 1.1, double tol=1.e-7, long iterMAX= 50);
    // destructor
    virtual ~ FD_vanilla(){};
    
    // methods
    
    double valuation(double S, double t) const;
    double delta(double S, double t) const;
    double gamma(double S, double t) const;
    double theta(double S, double t) const;
    
    int computationMethod() const   {return pdeKernal->computationMethod();};
    double computationAlpha() const {return pdeKernal-> computationAlpha();};
    bool computationSuccess() const {return pdeKernal-> computationSuccess();};
    
    
    HeatPDE* PDE_kernal() const{return new HeatPDE(*pdeKernal); };
    long logPriceDiscret() const {return N;};
    long expirationDiscret() const {return M;};
    option optionInfo() const {return opt;};
};


#endif /* FD_vanilla_hpp */
