#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h> 
//[[Rcpp::depends(RcppParallel)]]
#include <list>
#include <iterator>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace Rcpp; 
using namespace RcppParallel;
#include "linkPlink.h"
#include "glmmaiUpdate_plink.h"


//' Variance component update per iteration
//'
//' @title Update variance components tau in GLMM once
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param W 			numeric vector.
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param tol         	Numeric. Tolerance for minimum of variance components.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A List with the following elements:
//' * tau
//' * cov = XSiX_inv
//' * XSiX_inv_SiXt
//' * alpha
//' * eta
//' * Sigma_iX
// [[Rcpp::export]]
List glmmaiUpdate_plink(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::sp_mat& sparsekin, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){
	// Function warning("warning");
	List re(7);
	int n = n_nomissing;
	int ncol_x = X.n_cols;
	arma::fmat Sigma_iX(n,ncol_x);

	int q2 = arma::sum(fixtau<1);
	arma::uvec idxtau = arma::find(fixtau<1);

	re = getScore_plink(n_nomissing, q, Y, X, sparsekin, W, tau, fixtau, TRUE, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);

	// update tau: 2/2 (1/2 in getScore_plink)
	arma::fvec dtau = re["dtau"];	
	Rcout << "dtau estimates: " << dtau << std::endl;
	arma::fvec tau_old = tau, tau_new = tau;

	if(q2>0){
		// tau_new[idxtau] = tau_old[idxtau] + dtau[idxtau];
		for(int iter=0; iter < q2; iter++){
			tau_new( idxtau(iter) ) = tau_old( idxtau(iter) ) + dtau( idxtau(iter) );
		}	
		float step = 1.0;
	  	while ( arma::any(tau_new < 0) ){ 
			step = step*0.5; 
			// tau_new[idxtau] = tau_old[idxtau] + step*dtau[idxtau];
			for(int iter=0; iter < q2; iter++){
				tau_new( idxtau(iter) ) = tau_old( idxtau(iter) ) + step * dtau( idxtau(iter) );
				if(tau_new( idxtau(iter) )<tol && tau_old( idxtau(iter) )<tol){
					tau_new( idxtau(iter) ) = 0;
				}
			}	
	 	} // end while	
	 	tau_new.elem( arma::find(tau_new < tol) ).zeros(); 
	} // end if q2

 	if(tau_new(0)<tol && fixtau(0)<1){
		tau_new(0) = tau_old(0);
	}
    
 	List out = List::create(Named("tau") = tau_new, Named("cov") = re["cov"], Named("XSiX_inv_SiXt")=re["XSiX_inv_SiXt"], Named("alpha") = re["alpha"], Named("eta") = re["eta"]);
	return(out);


}

//' Variance component update per iteration
//' without sparse kinship matrix
//'
//' @title Update variance components tau in GLMM once
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param W 			numeric vector.
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param tol         	Numeric. Tolerance for estimation of variance components.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A List with the following elements:
//' * tau
//' * cov = XSiX_inv
//' * XSiX_inv_SiXt
//' * alpha
//' * eta
//' * Sigma_iX
// [[Rcpp::export]]
List glmmaiUpdate_plink_nok(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){
	// Function warning("warning");
	List re(7);
	int n = n_nomissing;
	int ncol_x = X.n_cols;
	arma::fmat Sigma_iX(n,ncol_x);

	int q2 = arma::sum(fixtau<1);
	arma::uvec idxtau = arma::find(fixtau<1);

	re = getScore_plink_nok(n_nomissing, q, Y, X, W, tau, fixtau, TRUE, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);
	
	// update tau: 2/2 (1/2 in getScore_plink)
	arma::fvec dtau = re["dtau"];	
	Rcout << "dtau estimates: " << dtau << std::endl;
	arma::fvec tau_old = tau, tau_new = tau;

	if(q2>0){
		for(int iter=0; iter < q2; iter++){
			tau_new( idxtau(iter) ) = tau_old( idxtau(iter) ) + dtau( idxtau(iter) );
		}
		float step = 1.0;
	  	while ( arma::any(tau_new < 0) ){ 
			step = step*0.5; 
			// tau_new[idxtau] = tau_old[idxtau] + step*dtau[idxtau];
			for(int iter=0; iter < q2; iter++){
				tau_new( idxtau(iter) ) = tau_old( idxtau(iter) ) + step * dtau( idxtau(iter) );
				if(tau_new( idxtau(iter) )<tol && tau_old( idxtau(iter) )<tol){
					tau_new( idxtau(iter) ) = 0;
				}
			}	
	 	} // end while	
	 	tau_new.elem( arma::find(tau_new < tol) ).zeros(); 		
	} // end if q2
	if(tau_new(0)<tol && fixtau(0)==0){
		tau_new(0) = tau_old(0);
	}

 	List out = List::create(Named("tau") = tau_new, Named("cov") = re["cov"], Named("XSiX_inv_SiXt")=re["XSiX_inv_SiXt"], Named("alpha") = re["alpha"], Named("eta") = re["eta"]);
	return(out);


}


//' Variance component update per iteration (LOCO)
//'
//' @title Update variance components tau in GLMM once
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param W 			numeric vector.
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param tol         	Numeric. Tolerance for estimation of variance components.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A List with the following elements:
//' * tau
//' * cov = XSiX_inv
//' * XSiX_inv_SiXt
//' * alpha
//' * eta
//' * Sigma_iX
// [[Rcpp::export]]
List glmmaiUpdate_plink_LOCO(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::sp_mat& sparsekin, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){
	// Function warning("warning");
	List re(7);
	int n = n_nomissing;
	int ncol_x = X.n_cols;	
	arma::fmat Sigma_iX(n,ncol_x);

	int q2 = arma::sum(fixtau<1);
	arma::uvec idxtau = arma::find(fixtau<1);
		
	re = getScore_plink_LOCO(n_nomissing, q, Y, X, sparsekin, W, tau, fixtau, TRUE, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);
	
	// update tau: 2/2 (1/2 in getScore_plink)
	arma::fvec dtau = re["dtau"];	
	Rcout << "dtau estimates: " << dtau << std::endl;
	arma::fvec tau_old = tau, tau_new = tau;
	// tau_new = tau_old + dtau;
	if(q2>0){
		for(int iter=0; iter < q2; iter++){
			tau_new( idxtau(iter) ) = tau_old( idxtau(iter) ) + dtau( idxtau(iter) );
		}		
	} // end if q2
	if(tau_new(0)<tol && fixtau(0)==0){
		tau_new(0) = tau_old(0);
	}

	// // set tau values less than tol to 0
	// tau_new.elem( arma::find(tau_new < tol) ).zeros(); 
  	
  	float step = 1.0;
  	// adjust negative values to 0
  	while ( arma::any(tau_new < 0.0) ){ 
		step = step*0.5;
		tau_new = tau_old + step * dtau; 
		// tau_new.elem( arma::find(zero_vec % (tau_new < tol)) ).zeros(); 
		// tau_new.elem( arma::find((tau_new < tol) and (tau_old < tol)) ).zeros(); 
		tau_new.elem( arma::find(tau_new < tol) ).zeros(); 

 	}
 	// set tau values less than tol to 0
 	tau_new.elem( arma::find(tau_new < tol) ).zeros();
 // 	if(tau_new(0)<tol && fixtau(0)==0){
	// 	tau_new(0) = tau_old(0);
	// }


 	List out = List::create(Named("tau") = tau_new, Named("cov") = re["cov"], Named("XSiX_inv_SiXt")=re["XSiX_inv_SiXt"], Named("alpha") = re["alpha"], Named("eta") = re["eta"]);
	return(out);


}
