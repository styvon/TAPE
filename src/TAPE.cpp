// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "TAPE.hpp"
#include "approxfun.hpp"
#include <math.h>   


TAPE::TAPEClass::TAPEClass(arma::mat t_cumul,
                         arma::vec t_mresid,
                         int t_N,
                         double t_SPA_Cutoff,
                         double t_var_ratio)
{
  m_mresid = t_mresid;
  m_varResid = var(m_mresid);
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_var_ratio = t_var_ratio;
  
  m_K_0_emp.setApproxFun(t_cumul.col(0), t_cumul.col(1));
  m_K_1_emp.setApproxFun(t_cumul.col(0), t_cumul.col(2));
  m_K_2_emp.setApproxFun(t_cumul.col(0), t_cumul.col(3));
  

}

double TAPE::TAPEClass::K_0(double t, 
                      int N0, 
                      double adjG0, 
                      arma::vec adjG1)        // adjusted Genotype 
{
  double t_adjG0 = t * adjG0;
  arma::vec t_adjG1 = t * adjG1;
  double out = N0 * m_K_0_emp.getValue(t_adjG0) + arma::sum(m_K_0_emp.getVector(t_adjG1));
  return out;
}

double TAPE::TAPEClass::K_1(double t,
                      int N0, 
                      double adjG0, 
                      arma::vec adjG1,        // adjusted Genotype
                      double q2)
{
  double t_adjG0 = t * adjG0;
  arma::vec t_adjG1 = t * adjG1;
  double out = N0 * adjG0 * m_K_1_emp.getValue(t_adjG0) + arma::sum(adjG1 % m_K_1_emp.getVector(t_adjG1)) - q2;
  return out;
}

double TAPE::TAPEClass::K_2(double t, 
                      int N0, 
                      double adjG0, 
                      arma::vec adjG1)       // adjusted Genotype
{
  double t_adjG0 = t * adjG0;
  arma::vec t_adjG1 = t * adjG1;
  double out = N0 * pow(adjG0, 2) * m_K_2_emp.getValue(t_adjG0) + arma::sum(pow(adjG1, 2) % m_K_2_emp.getVector(t_adjG1));
  return out;
}

Rcpp::List TAPE::TAPEClass::fastgetroot_K1(double t_initX,
                                     int N0, 
                                     double adjG0, 
                                     arma::vec adjG1,        // adjusted Genotype
                                     double q2)
{
  double x = t_initX, oldX;
  double K1 = 0, K2 = 0, oldK1;
  double diffX = arma::datum::inf, oldDiffX;
  bool converge = true;
  double tol = 0.001;
  int maxiter = 100;
  int iter = 0;
  
  for(iter = 0; iter < maxiter; iter ++){
    
    oldX = x;
    oldDiffX = diffX;
    oldK1 = K1;
    
    K1 = K_1(x, N0, adjG0, adjG1, q2);
    K2 = K_2(x, N0, adjG0, adjG1);
    
    diffX = -1 * K1 / K2;
    
    if(!std::isfinite(K1)){
      // checked it on 07/05:
      // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
      // then we can set the one sided p value as 0 (instead of setting converge = F)
      x = arma::datum::inf;
      K2 = 0;
      break;
    }
    
    if(arma::sign(K1) != arma::sign(oldK1)){
      while(std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }
    
    if(std::abs(diffX) < tol) break;
    
    x = oldX + diffX;
  }
  
  if(iter == maxiter) 
    converge = false;
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                        Rcpp::Named("iter") = iter,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K2") = K2);
  return yList;
}

double TAPE::TAPEClass::GetProb_SPA(double adjG0, 
                   arma::vec adjG1,
                   int N0, 
                   double q2, 
                   bool lowerTail)
{
  double initX = 0;
  
  // The following initial values are validated on 03/25/2021
  if(q2 > 0) initX = 3;
  if(q2 <= 0) initX = -3;
  
  Rcpp::List rootList = fastgetroot_K1(initX, N0, adjG0, adjG1, q2);
  double zeta = rootList["root"];
  
  double k1 = K_0(zeta,  N0, adjG0, adjG1);
  double k2 = K_2(zeta,  N0, adjG0, adjG1);
  double temp1 = zeta * q2 - k1;
  
  double w = arma::sign(zeta) * sqrt(2 * temp1);
  double v = zeta * sqrt(k2);
  
  double pval = arma::normcdf(arma::sign(lowerTail-0.5) * (w + 1/w * log(v/w)));
  
  return pval;
}

void TAPE::TAPEClass::getMarkerPval(arma::vec t_GVec, 
                                    double& t_Beta,
                                    double& t_seBeta, 
                                    double& t_pval, 
                                    double t_altFreq, 
                                    double& t_zScore)
{
  double S = sum(t_GVec % m_mresid);
  arma::vec adjGVec = t_GVec - 2 * t_altFreq;
  arma::vec adjGVec2 = pow(adjGVec, 2);
  double VarS = m_varResid * sum(adjGVec2) * m_var_ratio;
  
  t_zScore = S / sqrt(VarS);
  double pvalNorm = 2 * arma::normcdf(-1*t_zScore);
  double pval = pvalNorm;
  
  if(std::abs(t_zScore) > m_SPA_Cutoff){
    arma::uvec N1set = arma::find(t_GVec!=0);  // position of non-zero genotypes
    int N0 = m_N - N1set.size();
    
    arma::vec adjGVecNorm = adjGVec / sqrt(VarS); // normalized genotype (such that sd=1)
    
    arma::vec adjG1 = adjGVecNorm.elem(N1set);
    double adjG0 = -2 * t_altFreq / sqrt(VarS);  // all subjects with g=0 share the same normlized genotype, this is to reduce computation time
    
    double pval1 = GetProb_SPA(adjG0, adjG1, N0, std::abs(t_zScore), false);
    double pval2 = GetProb_SPA(adjG0, adjG1, N0, -1*std::abs(t_zScore), true);
    double pval = pval1 + pval2;
    
  } // end if(std::abs(t_zScore) > m_SPA_Cutoff)
  
  t_pval = pval;
  t_Beta = S / VarS;
  t_seBeta = t_Beta / t_zScore;
}
