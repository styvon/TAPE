#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
#include <list>
#include <iterator>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string> // for string class



using namespace Rcpp; 
using namespace std;
using namespace arma;

// incomplete cholesky factorization
// [[Rcpp::export]]
arma::fmat icc(arma::fmat A){
  int N = A.n_cols ;
  arma::fmat temp = A;
  for(int k = 0; k < N; k++){
    temp(k,k) = sqrt(temp(k,k));
    for(int i = k + 1; i < N; i++){
      if(temp(i,k) != 0){
        temp(i,k) = temp(i,k)/temp(k,k);
      }
    }
    for(int j = k + 1; j < N; j++){
      for(int i= j; i < N; i++){
        if(temp(i,j) != 0){
          temp(i,j) = temp(i,j) - temp(i,k)*temp(j,k);
        }
      }
    }
  }
  
  for(int i = 0; i<N; i++){
    for(int j = i+1; j<N; j++){
      temp(i,j) = 0;
    }
  }
  
  return temp;
}


void set_seed(unsigned int seed) {
	// R func for set seed
	Rcpp::Environment base_env("package:base");
  	Rcpp::Function set_seed_r = base_env["set.seed"];
  	set_seed_r(seed);  
}


NumericVector nb(int n) { 
	// R func for generate binomial rv
  	return(rbinom(n,1,0.5));
}


float calCV(arma::fvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = arma::mean(xVec);
  float vecSd = arma::stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}