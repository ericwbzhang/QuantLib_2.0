//
//  lm.hpp
//  QuantLib_2.0
//
//  Created by EricZ on 2/13/16.
//  Copyright © 2016 EricZ. All rights reserved.
//

#ifndef lm_hpp
#define lm_hpp

#include <iostream>
#include <stdio.h>
#include "algebra.h"
#include <Boost/math/distributions.hpp>
#include <random>
#include <Math/IFunction.h>
#include <Math/BrentMinimizer1D.h>


// This module contains several linear regression models. They are all based on least square method.



class lm{
    
protected:
    double R2;
    double AdjR2;
    double res_variance;
    long n, p;
    double df;
    
    Eigen::VectorXd coef;
    Eigen::VectorXd res;
    Eigen::MatrixXd t_stats;
    Eigen::MatrixXd F_stats;
    Eigen::MatrixXd coef_cov;
    
    virtual void computeStats(){
        // AdjR2
        AdjR2=1- (1-R2)*(n-1)/double(n-df);
        
        // t stats
        t_stats.resize(p,2);
        boost::math::students_t_distribution<double > t_dist(fmax(long (n- df),1));
        for (long i=0; i<p; i++){
            t_stats(i,0)= coef(i)/ sqrt(coef_cov(i,i));
            t_stats(i,1)= 1- cdf(t_dist, t_stats(i,0));
        }
        //F stats
        F_stats.resize(1,2);
        F_stats(0,0)= R2/(1-R2)* double(n- df)/ (df-1);
        boost::math::fisher_f_distribution<double> f_dist(fmax(long (df-1), 1), long (n-df));
        F_stats(0,1)= 1-cdf(f_dist, fmax(0.0,F_stats(0,0)));
        
    };
    
    std::string coutAdding;
public:
    lm(){};
    virtual ~lm(){};
    
    virtual double Rsquare(){return R2;};
    virtual double resVariance(){return res_variance;};
    virtual Eigen::MatrixXd coefCovariance(){return coef_cov; };
    virtual Eigen::MatrixXd tStats(){ return t_stats; };
    virtual Eigen::MatrixXd FStats(){ return F_stats;} ;
    virtual Eigen::VectorXd coeff(){return coef;};
    virtual Eigen::VectorXd residual(){return res;};
    virtual long obs_count(){return n;};
    virtual long regressor_count(){return p;};
    virtual long degreeFreedom(){return df;};
    
    virtual Eigen::VectorXd prediction(const Eigen::MatrixXd & X){
        return X* coef;
    };
    virtual Eigen::VectorXd predictionError( const Eigen::MatrixXd &X, const Eigen::VectorXd &y){
        return y- X*coef;
    };
    virtual double predictionErrorVar(const Eigen::MatrixXd & X, const Eigen::VectorXd & y){
        Eigen::VectorXd error(X*coef-y);
        return (error.squaredNorm()-error.size()* pow(error.mean(),2.0))/ (error.size()-1);
    };
    virtual double predictionErrorStdiv(const Eigen::MatrixXd &X, const Eigen::VectorXd &y){
        return sqrt(this-> predictionErrorVar(X, y));
    }
    
    
    
    friend std::ostream & operator<<( std::ostream &output,const lm & reg  ){
        
        Eigen::MatrixXd tmp( reg.t_stats.rows(),reg.t_stats.cols()+1);
        tmp.col(0)= reg.coef;
        tmp.block(0,1,tmp.rows(), tmp.cols()-1)= reg.t_stats;
        
        output << "Linear Model:  \n*********************\n" << reg.coutAdding<<reg.n<< " observations, " << reg.p<< " regressors (including intercept). Degree Freedom: "<< reg.df<<" \nR square: "<< reg.R2<< "\t Adjusted R square: "<< reg.AdjR2<<std::endl<< "F Stats: \n"<< reg.F_stats<< std::endl<<"\nCoefficients\nEst\t\t\tt_Stats\t\t\tpValue\n"<< tmp<< std::endl<< std::endl<< "Residual stdiv: "<<  sqrt(reg.res_variance)<< std::endl<< std::endl;
        return output;
    };
};


// lm_OLS
class lm_OLS: public lm {
    
public:
    lm_OLS(){};
    virtual ~lm_OLS(){};
    lm_OLS(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, bool addConst=false){
        
        // adjust X to accomdate addConst
        Eigen::MatrixXd A(X);
        if(addConst) {
            long nrow= A.rows();
            long ncol= A.cols();
            Eigen::VectorXd ones;  ones.setOnes(nrow);
            A.resize(nrow, ncol+1);
            A.col(0)= ones;
            A.block(0,1, nrow, ncol)= X;
            
        }
        n= A.rows();
        p= A.cols();
        df=p;
        // estimate the coefs
        SVD_Jb svd(A, Eigen::ComputeThinU, Eigen::ComputeThinV);
        coef= svd.solve(y);
        // compute residuals and residual variance
        res= y- A*coef;
        res_variance= res.squaredNorm()/ (n- df);
        // estimate coef covariance
        Eigen::VectorXd S_diag= svd.matrixS().diagonal();
        long m= S_diag.size();
        for (long i=0; i<m; i++) S_diag(i)= 1.0/pow(S_diag(i),2.0);
        coef_cov= svd.matrixV()* S_diag.asDiagonal() *svd.matrixV().transpose() * res_variance;
        // compute R square
        double y_variation= y.squaredNorm()-n* pow(y.mean(), 2.0);
        R2= 1- res.squaredNorm()/ y_variation;
        
        this->computeStats();
        
        coutAdding="OLS Regression\n\n";
    };
};


// lm_feasibleGLS
class lm_feasibleGLS: public lm{
    
protected:
    bool success; //whether convergence or not
    long iterations;
    
public:
    lm_feasibleGLS(){};
    virtual ~lm_feasibleGLS(){};
    lm_feasibleGLS(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, const Eigen::MatrixXd &Omega= Eigen::MatrixXd(), double tol=1e-5, long MaxIteration=50, bool addConst=false){
        // adjust X to accomdate addConst
        Eigen::MatrixXd A(X);
        if(addConst) {
            long nrow= A.rows();
            long ncol= A.cols();
            Eigen::VectorXd ones;  ones.setOnes(nrow);
            A.resize(nrow, ncol+1);
            A.col(0)= ones;
            A.block(0,1, nrow, ncol)= X;
            
        }
        n= A.rows();
        p= A.cols();
        df=p;
        // estimate the coefs
        Eigen::MatrixXd W; W.setIdentity(n,n);
        if(Omega.rows()== n && Omega== Omega.transpose() && Omega.determinant()!=0) {
            W= Omega;
        }
        W= W.inverse();
        coef= (A.transpose()* W* A).inverse() *A.transpose() *W *y;
        Eigen::VectorXd res_old= y- A*coef;
        
        long count=1;
        double error= 1000; //error tracks the  normal of (res_new- res_old)
        success= false;
        while (error> tol && count< MaxIteration){
            Eigen::VectorXd tmp(res_old);
            for (long i=0; i<n ; i++) tmp(i)= 1/pow(tmp(i),2.0);
            W= tmp.asDiagonal();
            coef= (A.transpose()* W* A).inverse() *A.transpose() *W *y;
            res= y- A*coef;
            error= (res- res_old).norm();
            count++;
            res_old= res;
        }
        if (error< tol) success= true;
        iterations= count;
        // compute residual variance
        res_variance= res.squaredNorm()/ (n- df);
        // estimate coef covariance
        coef_cov=  (A.transpose() * W * A).inverse() ;
        // compute R square
        double y_variation= y.squaredNorm()-n* pow(y.mean(), 2.0);
        R2= 1- res.squaredNorm()/ y_variation;
        
        this->computeStats();
        
        std::ostringstream ss; ss <<"feasible GLS Regression. Success= "<< success << " with "<<iterations<<" iterations. \n\n";
        coutAdding= ss.str();
        
    };
};


// Weighted LS

class lm_WLS: public lm{
    
public:
    lm_WLS(){};
    virtual ~lm_WLS(){};
    lm_WLS(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, const Eigen::MatrixXd &W){
        // we require W to be symmetric
        
        // adjust X to accomdate addConst
        Eigen::MatrixXd A(X);
        n= A.rows();
        p= A.cols();
        df=p;
        // estimate the coefs
        Eigen::MatrixXd tmp= (A.transpose()* W* A).inverse();
        coef= tmp* A.transpose()* W* y;
        // compute residuals and residual variance
        res= y- A*coef;
        res_variance= res.squaredNorm()/ (n- df);
        // estimate coef covariance
        coef_cov= tmp* A.transpose()* W* W.transpose()* A* tmp * res_variance;
        // compute R square
        double y_variation= y.squaredNorm()-n* pow(y.mean(), 2.0);
        R2= 1- res.squaredNorm()/ y_variation;
        
        this->computeStats();
        
        coutAdding="Weighted LS Regression\n\n";
    };
    
    
};

// Ridge Regression

class lm_Ridge: public lm{
    
protected:
    double lambda;
    long k;
    
    
    class CV_error: public ROOT::Math::IBaseFunctionOneDim{
    private:
        std::vector<Eigen::MatrixXd> TestX;
        std::vector<Eigen::VectorXd> TestY;
        std::vector<Eigen::MatrixXd> TrainX;
        std::vector<Eigen::VectorXd> TrainY;
        
        long N;
        
        
        void ithFold( const Eigen::MatrixXd &A, const Eigen::VectorXd &y, Eigen::MatrixXd &trainX, Eigen::VectorXd &trainY, Eigen::MatrixXd &testX, Eigen::VectorXd &testY, long k, long i ){
            // jthFold generates the jth training and testing sample of the k-fold validation
            long n= A.rows();
            long p= A.cols();
            long test_size= n/k;
            i=i-1;
            // construct traning sample and test sample
            
            long size; // size defines how many rows to include in test sample
            long test_start; //  the row index of the starting test sample
            long test_end; // the row index of the ending test sample
            
            if(i==k-1){
                size= n-(k-1)*test_size;
                test_start= i*test_size;
                
                testX= A.block(test_start, 0, size, p);
                testY= y.block(test_start, 0, size, 1);
                trainX= A.block(0,0, n-size, p);
                trainY= y.block(0,0, n-size, 1);
                
            }else {
                size= test_size;
                test_start= i* test_size;
                test_end= test_start+ size-1;
                
                testX= A.block(test_start, 0, size, p);
                testY= y.block(test_start, 0, size, 1);
                trainX.resize(n-size, p);
                trainY.resize(n-size);
                trainX.block(0,0, test_start, p)= A.block(0,0, test_start, p);
                trainX.block(test_start,0, n-size- test_start,p)= A.block(test_end+1, 0, n- size- test_start, p);
                trainY.block(0,0, test_start,1)= y.block(0,0, test_start,1);
                trainY.block(test_start,0, n-size- test_start,1)= y.block(test_end+1,0, n-size- test_start,1);
            }
            
        };
        
        double DoEval(double lam) const{
            double error=0;
            
            long n= TrainX.size();
            for (long i=0; i<n; i++){
                // train the regression model with lambda= lam
                lm_Ridge lmR(TrainX[i], TrainY[i], lam);
                // evaluate and record prediction error
                Eigen::VectorXd tmp= lmR.predictionError(TestX[i], TestY[i]);
                
                error+= tmp.squaredNorm();
            }
            
            error/= N;
            return error;
        };
        
    public:
        CV_error(){};
        CV_error(const Eigen::MatrixXd &A, const Eigen::VectorXd &y, long k){
            N= A.rows();
            for (long i=1 ;i <=k ; i++){
                Eigen::MatrixXd trainX, testX;
                Eigen::VectorXd trainY, testY;
                this->ithFold(A, y, trainX, trainY, testX, testY, k, i);
                TrainX.push_back(trainX);
                TrainY.push_back(trainY);
                TestX.push_back(testX);
                TestY.push_back(testY);
            }
        };
        virtual ~CV_error(){};
        
        ROOT::Math::IBaseFunctionOneDim * Clone() const { return new CV_error(*this);};
        
    };
    
    bool optiLambda(const Eigen::MatrixXd &A, const Eigen::VectorXd &y, double lamMIN, double lamMAX, double tol1, double tol2, long iterMAX, double & lam, std::stringstream &log ){
        // optiLambda wraps the ROOT::Math::BrentMinimizer1D minimizer. It returns a flag whether the minimizer converges.
        // tol1 controls the abs error; tol2 controls the relative error
        
        CV_error CVErr(A, y, k);
        ROOT::Math::BrentMinimizer1D bm;
        bm.SetFunction(CVErr, lamMIN, lamMAX);
        bool res= bm.Minimize(iterMAX, tol1, tol2);
        
        lam= bm.XMinimum();
        
        log<< "Convergence: " << res<< std::endl;
        return res;
    };
    
public:
    lm_Ridge(){};
    virtual ~lm_Ridge(){};
    lm_Ridge(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, double lam, bool addConst=false){
        lambda= lam;
        k=0;
        
        // adjust X to accomdate addConst
        Eigen::MatrixXd A(X);
        if(addConst) {
            long nrow= A.rows();
            long ncol= A.cols();
            Eigen::VectorXd ones;  ones.setOnes(nrow);
            A.resize(nrow, ncol+1);
            A.col(0)= ones;
            A.block(0,1, nrow, ncol)= X;
            
        }
        n= A.rows();
        p= A.cols();
        
        
        // estimate the coefs
        Eigen::MatrixXd tmp; tmp.setIdentity(p,p);
        tmp=(A.transpose()* A+ lambda*tmp).inverse();
        coef= tmp* A.transpose()* y;
        // define df
        df= (A* tmp* A.transpose()).trace();
        
        // compute residuals and residual variance
        res= y- A*coef;
        res_variance= res.squaredNorm()/ (n- df);
        // estimate coef covariance
        coef_cov= tmp* A.transpose()* A* tmp* res_variance;
        // compute R square
        double y_variation= y.squaredNorm()-n* pow(y.mean(), 2.0);
        R2= 1- res.squaredNorm()/ y_variation;
        
        this->computeStats();
        
        std::ostringstream ss; ss <<"Ridge Regression. Lambda= "<< lambda << ". Degree Freedom= "<<df<<"\nNote: For Ridge regression, both F-Stats and t-Stats and R^2 are the rough. (Since the df is continuous) \n\n";
        coutAdding= ss.str();
        
    };
    
    
    
    lm_Ridge(const Eigen::MatrixXd &X, const Eigen::VectorXd &y,
             std::vector<double> lamRange={0,100}, double tol_multiplier= 1e-3, double tol_error= 1e-4, long iterMAX= 100,
             long K=10, bool addConst= false){
        // lamRange= {lamMIN, lamMAX}. lamMax and lamMin set the range of lambda. k sets how many folds we do in cross validation.
        // tol_multiple sets the stopping condition of optimal lambda search. We stop is the interval length are smaller than tol_multiple*(lamMAX- lamMIN).
        // tol_error controls the change of error. If error updates less than tol_error we stop evolve lambda.
        // iterMAX controls the max of iteration times.
        
        
        k=K;
        // adjust X to accomdate addConst
        Eigen::MatrixXd A(X);
        if(addConst) {
            long nrow= A.rows();
            long ncol= A.cols();
            Eigen::VectorXd ones;  ones.setOnes(nrow);
            A.resize(nrow, ncol+1);
            A.col(0)= ones;
            A.block(0,1, nrow, ncol)= X;
            
        }
        n= A.rows();
        p= A.cols();
        
        // decides optimal lambda.
        
        // randomly permutate A and y
        //first randomly permutate A and y
        Eigen::PermutationMatrix<Eigen::Dynamic> perm;
        perm.setIdentity(n);
        std::shuffle(perm.indices().data(), perm.indices().data()+ perm.indices().size(), std::default_random_engine(int(time(0))));
        Eigen::MatrixXd permA= perm*A;
        Eigen::VectorXd permy= perm*y;
        // define CVerror
        std::stringstream log;
        this->optiLambda(permA, permy, lamRange[0], lamRange[1], tol_error, tol_multiplier*(lamRange[1]-lamRange[0]), iterMAX, lambda, log);
        
        // estimate the coefs
        Eigen::MatrixXd tmp; tmp.setIdentity(p,p);
        tmp=(A.transpose()* A+ lambda*tmp).inverse();
        coef= tmp* A.transpose()* y;
        // define df
        df= (A* tmp* A.transpose()).trace();
        
        // compute residuals and residual variance
        res= y- A*coef;
        res_variance= res.squaredNorm()/ (n- df);
        // estimate coef covariance
        coef_cov= tmp* A.transpose()* A* tmp* res_variance;
        // compute R square
        double y_variation= y.squaredNorm()-n* pow(y.mean(), 2.0);
        R2= 1- res.squaredNorm()/ y_variation;
        
        this->computeStats();
        
        std::ostringstream ss; ss <<"Ridge Regression. Lambda= "<< lambda << ". Degree Freedom= "<<df<<"\nNote: For Ridge regression, both F-Stats and t-Stats and R^2 are the rough. (Since the df is continuous) \n\n";
        coutAdding= ss.str()+ log.str()+"\n\n";
        
    };
    
    
};


/*
 ************ test *************
 
 int main(){
 
 //    std::srand(int(time(0)));
 //    Eigen::MatrixXd A;
 //    A.setRandom(1000, 1)*1e4;
 //    // Eigen::MatrixXd C= A+ A.transpose();
 //    Eigen::MatrixXd B;
 //    B.setRandom(1000,1);B= 2*A+B;
 //    Eigen::MatrixXd C(1000,2); C<< A, A;
 //    C.col(0).setOnes();
 //
 //    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
 //    std::cout<< qr.solve(B)<<std::endl<<std::endl;
 //
 //    lm_Ridge lm(C,B);
 //    std::cout<< lm<<std::endl;
 
 return 0;
 }
 */
#endif /* lm_hpp */