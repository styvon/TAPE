#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::max


using namespace Rcpp; 
#include "pcgsolve.h"


//' Get trace using Hutchinson’s randomized trace estimator
//' @param n_nomissing  integer
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param densekin 	list of dense covariance matrices read from kmatfile.
//' @param Sigma 		matrix
//' @param Sigma_iX 	matrix
//' @param X 			matrix
//' @param W_vec 		vector of type arma::fvec.
//' @param XSiX_inv 	matrix, covariance
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
// [[Rcpp::export]]
arma::fvec getTrace(int n_nomissing, int q, arma::sp_mat& sparsekin, const std::vector<arma::fmat>& densekin, arma::fmat& Sigma, arma::fmat& Sigma_iX, arma::fmat& X, const arma::fvec& W_vec, arma::fmat& XSiX_inv,  int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace){
  	set_seed(200);
  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = n_nomissing;
  	arma::fmat temp_mat(nrun_trace, q+1);
  	arma::fvec temp_vec(nrun_trace);
    arma::vec temp_vec_double(Nnomissing);
  	temp_mat.zeros();

    arma::fvec Sigma_iu;
    arma::fvec Pu;
    arma::fmat Au_mat(Nnomissing, q+1);
    arma::fvec u_vec;
    NumericVector u_vec0;

    int nrun_trace_start = 0;
    int nrun_trace_end = nrun_trace;
    arma::fvec trace_cv(q+1);
    trace_cv.fill(cutoff_trace + 0.1);

    // while((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
	while( arma::any(trace_cv > cutoff_trace) ){

    	for(int i = nrun_trace_start; i < nrun_trace_end; i++){

    		u_vec0 = nb(Nnomissing);
    		u_vec = as<arma::fvec>(u_vec0);
    		u_vec = u_vec*2 - 1;
 
    		// Sigma_iu = pcgsolve(Sigma, u_vec, "Jacobi", tol_pcg, maxiter_pcg);
        Sigma_iu = arma::solve(Sigma, u_vec, arma::solve_opts::likely_sympd);


    		Pu = Sigma_iu - Sigma_iX * (XSiX_inv *  (Sigma_iXt * u_vec));
        Au_mat.col(0) = u_vec;
        temp_mat(i,0) = dot(Au_mat.col(0), Pu);

        // conversion for ops with sp_mat
        temp_vec_double = 0.0+sparsekin * arma::conv_to<arma::vec>::from(u_vec);
        Au_mat.col(1) = arma::conv_to<arma::fvec>::from(temp_vec_double);
        // Au_mat.col(1) = sparsekin * u_vec;
        temp_mat(i,1) = dot(Au_mat.col(1), Pu);

    		for(int j=2; j<q+1;j++){
    			Au_mat.col(j) = densekin[j-2] * u_vec;
    			temp_mat(i,j) = dot(Au_mat.col(j), Pu);
    		} // end for j in 2:q
    		
  		} // end for i


  		// update trace cv vector
  		for(int i=0; i<q+1; i++){
  			temp_vec = temp_mat.col(i);
  			trace_cv(i) = calCV( temp_vec );
  		} 

  		// if not converge, increase nrun_trace and rerun
  		if( arma::any(trace_cv > cutoff_trace) ){
  		// if((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
  			nrun_trace_start = nrun_trace_end;
  			nrun_trace_end = nrun_trace_end + 10;
  			temp_mat.reshape(nrun_trace_end,q+1);
  			//std::cout << "arma::mean(temp_mat0): " << arma::mean(temp_mat0) << std::endl;	
  			Rcout << "CV for trace random estimator using "<< nrun_trace_start << " runs is " << trace_cv <<  "(> " << cutoff_trace << std::endl;
  			Rcout << "try " << nrun_trace_end << "runs" << std::endl;
      } // end if arma::any(trace_cv > cutoff_trace)

    } // end while  arma::any(trace_cv > cutoff_trace)
    Au_mat.clear();
    Pu.clear();
    Sigma_iu.clear();
    u_vec.clear();
    temp_vec.clear();

  	arma::fvec traVec(q+1);
  	for(int i=0; i<q+1; i++){
  		traVec(i) = arma::mean(temp_mat.col(i));
  	}
    temp_mat.clear();
	
  	return(traVec);
}