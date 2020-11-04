#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
#include <list>
#include <iterator>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string> // for string class

#include "pcgUtil.h"

using namespace Rcpp; 
using namespace std;
using namespace arma;



//' Preconditioned conjugate gradient method (from cPCG)
//'
//' Preconditioned conjugate gradient method for solving system of linear equations Ax = b,
//' where A is symmetric and positive definite.
//'
//' @title Solve for x in Ax = b using preconditioned conjugate gradient method.
//' @param A                matrix, symmetric and positive definite.
//' @param b                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
// [[Rcpp::export]]
arma::fvec pcgsolve(arma::fmat & A, arma::fvec & b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000) {

  // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;
  
  // get preconditioner M
  arma::fmat M;
  if (preconditioner == "Jacobi"){
    M = arma::diagmat(A);
  } else if(preconditioner == "SSOR"){
    arma::fmat D = arma::diagmat(A);
    arma::fmat L = arma::trimatl(A);
    M = (D+L) * D.i() * (D+L).t();
  } else if(preconditioner == "ICC"){
    M = icc(A);
  }

  
  // initiate solution x as zeros
  arma::fvec x(C) ;
  x.zeros() ; 
  
  arma::fvec oneVec(C);
  oneVec.ones() ;
  
  arma::fvec r = b - A * x;
  arma::fmat Minv = M.i();
  arma::fvec z = Minv * r;
  arma::fvec p = z;
  float rz_old = sum(r % z);
  float rz_new=1;
  // arma::fvec rz_ratio(1);

  
  arma::fvec Ap(R);
  float alpha, beta;
  // vector version of alpha
  // arma::fvec alphaVec(1);
  
  for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
    Ap = A * p;
    alpha = rz_old / sum(p % Ap);
    
    x += alpha * p;
    r -= alpha * Ap;
    z = Minv * r;
    rz_new = sum( z % r );
    beta = rz_new / rz_old; 
    
    p = z + beta * p;
    rz_old = rz_new;
    if (iter >= maxIter){
      Rcout << "pcg did not converge." << std::endl;
    }
  }
  
  return x;
  
} 
