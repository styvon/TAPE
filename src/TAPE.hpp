#ifndef TAPE_HPP
#define TAPE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "approxfun.hpp"

namespace TAPE{

class TAPEClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  approxfun::approxfunClass m_K_0_emp;
  approxfun::approxfunClass m_K_1_emp;
  approxfun::approxfunClass m_K_2_emp;
  
  arma::vec m_mresid;
  double m_varResid;
  arma::mat m_XinvXX, m_tX;
  int m_N;
  double m_pVal_covaAdj_Cutoff;
  double m_SPA_Cutoff;
  double m_var_ratio;
  
public:
  
  TAPEClass(arma::mat t_cumul,
            arma::vec t_mresid,
            int t_N,
            double t_SPA_Cutoff,
            double t_var_ratio);
  
  double K_0(double t, 
             int N0, 
             double adjG0, 
             arma::vec adjG1);
  
  double K_1(double t,
             int N0, 
             double adjG0, 
             arma::vec adjG1, 
             double q2);
  
  double K_2(double t, 
             int N0, 
             double adjG0, 
             arma::vec adjG1) ;
  
  Rcpp::List fastgetroot_K1(double t_initX,
                            int N0, 
                            double adjG0, 
                            arma::vec adjG1, 
                            double q2);
  
  double GetProb_SPA(double adjG0, 
                     arma::vec adjG1,
                     int N0, 
                     double q2, 
                     bool lowerTail);
  
  void getMarkerPval(arma::vec t_GVec, 
                     double& t_Beta, 
                     double& t_seBeta, 
                     double& t_pval, 
                     double t_altFreq, 
                     double& t_zScore);
  
};

}

#endif
