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
#include <string>


using namespace Rcpp; 
using namespace RcppParallel;
#include "glmmaiUpdate_plink.h"
#include "readMat.h"
#include "linkPlink.h"

//' Variance component estimation
//' Uses Plink file for GRM construction
//'
//' @title Iteratively estimate variance components tau in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param Y			column vector, working y
//' @param y			column vector, original y
//' @param X   			data matrix
//' @param W 			numeric vector.
//' @param alpha 		numeric vector.
//' @param maxiter 		integer, maximum iteration.
//' @param tol         	Numeric. Tolerance for convergence of variance components.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm. 
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm. 
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @param verbose		boolean. 
//' @param loco 		boolean. Indicator for leave-one-chromosome-out analysis.
//' @param startid_chr_LOCO integer vector. Indeces for starting point for 22 chromosomes.
//' @param endid_chr_LOCO 	integer vector. Indeces for ending point for 22 chromosomes.
//' @param tol_tau		Numeric. Threshold for minimum variance components value. Default is 1e-5.
//' @param tol_coef    Numeric. Tolerance for convergence of coefficients. Default is 0.1
//' @return A List.
//' * theta: variance components estimates
//' * coefficients
//' * linear_predictors
//' * fitted_values
//' * Y
//' * X
//' * residuals
//' * cov
//' * XSiX_inv_SiXt
//' * converged
//' * method
//' * LOCO
//' * result_LOCO (if LOCO=TRUE)
// [[Rcpp::export]]
List getTau_plink(int n_nomissing, arma::sp_mat& sparsekin, arma::fvec tau, arma::fvec fixtau, arma::fvec& Y, arma::fvec& y, arma::fmat& X, arma::fvec alpha, arma::fvec W, int maxiter, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace, bool verbose, bool loco, IntegerVector startid_chr_LOCO, IntegerVector endid_chr_LOCO, float tol_tau=1e-5, float tol_coef=0.1){

	int n =n_nomissing;
	// Kinships kinships(K, kmatfile);
	Function warning("warning");

	int q = 2; 
	
	arma::fvec W_new = W;
	// for iterative update
	arma::fvec alpha0, tau0;
	arma::fvec alpha_new = alpha, tau_new = tau;
	List fit(5);
	// arma::fmat cov(n,n);
	arma::fvec eta(n), mu(n), mu_eta(n);
	arma::fvec Y_new = Y; 

	int q2,i;
	arma::uvec idxtau = arma::find(fixtau<1);

	// for LOCO
	int startid, endid;
	LogicalVector startid_chr_LOCO_bool;
	LogicalVector endid_chr_LOCO_bool;
	bool startid_isna, endid_isna;
	List fit_LOCO(5);
	arma::fvec alpha0_LOCO, alpha_new_LOCO = alpha;
	arma::fvec eta_LOCO(n), mu_LOCO(n), mu_eta_LOCO(n);
	List results_LOCO(22);
	arma::fvec res_LOCO, Y_LOCO= Y;

	// ======== iterative update tau ===========
	for(i=0; i< maxiter; i++){
		if(verbose){
			Rcout << "------ Variance component estimation iteration " << i << " ------" << std::endl;
		}
		alpha0 = alpha_new;
		tau0 = tau_new;
		q2 = arma::sum(fixtau<1);

		// returns list with element "tau","alpha","eta","cov" (updated)
		fit = glmmaiUpdate_plink(n_nomissing, q, Y_new, X, W_new, sparsekin, tau_new, fixtau, tol_tau, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);	
		
		if(q2>0){
			arma::fvec temp1 = fit["tau"];
			tau_new = temp1;
			for(int j = 0; j<q+1; j++){
				if(fixtau(j)<1 && tau_new(j)< 1.01 * tol_tau){
					fixtau(j) = 1.0;
				}
			}
			
		} // end if q2

		arma::fvec temp2 = fit["alpha"];
		alpha_new = temp2;
		arma::fvec temp3 = fit["eta"];
		eta = temp3;


		if(verbose){
			Rcout << "Variance component estimates: " << tau_new << std::endl;
		}

		// for gaussian family
		mu = eta;
		mu_eta.fill(1.0);
		Y_new = eta + (y - mu)/mu_eta;
		// W is always ones for gaussian family, no need to update

		// check convergence
		arma::fvec check_vec_alpha = arma::abs(alpha_new - alpha0)/(arma::abs(alpha_new) + arma::abs(alpha0) + tol_coef);
		arma::fvec check_vec_tau = arma::abs(tau_new - tau0)/(arma::abs(tau_new) + arma::abs(tau0) + tol);

		if( (2*check_vec_tau.max() < tol) || ( 2*check_vec_alpha.max() < tol_coef) ){
			break;
		}

		float check_max = tau_new.max();
		if( check_max > 1.0/tol/tol) {
			warning("Large variance estimate observed in the iterations, model not converged...");
			i = maxiter;
			break;
		}
	} // end of iteration maxiter
	
	bool converged = (i < maxiter);
	arma::fvec res = y - mu; 

	arma::fmat Sigma_iX = pcgsolveMat_plink(sparsekin, W, tau_new, X, "Jacobi", tol_pcg, maxiter_pcg);
    
	if(!loco) {
		return List::create(Named("theta") = tau_new, Named("coefficients")=alpha_new, Named("linear_predictors")=eta, Named("fitted_values") = mu, Named("Y")=Y_new, Named("X")=X, Named("residuals")=res, Named("cov")=fit["cov"], Named("XSiX_inv_SiXt")=fit["XSiX_inv_SiXt"], Named("converged")=converged, Named("method") ="glmm_ai_multitau", Named("Sigma_iX")=Sigma_iX, Named("LOCO")=loco);
	}else{ 
		startid_chr_LOCO_bool = is_na(startid_chr_LOCO);
		endid_chr_LOCO_bool = is_na(endid_chr_LOCO);

		// ======== for each chromosome, iterative update tau (if LOCO) ===========
		for(int chr=1; chr<23; chr++){
			
			startid_isna = startid_chr_LOCO_bool[chr-1];			
			endid_isna = endid_chr_LOCO_bool[chr-1];
			if((!startid_isna) && (!endid_isna)){
				startid = startid_chr_LOCO[chr-1];
				endid = endid_chr_LOCO[chr-1];
				Rcout << "----- LOCO: chromosome " << chr << " -----" << std::endl;
				setLOCOid(startid, endid);
				for( i=0; i< maxiter; i++){
					if(verbose){
						Rcout << "------ Variance component estimation iteration " << i << " ------" << std::endl;
					}
					alpha0_LOCO = alpha_new_LOCO;

					fit_LOCO = glmmaiUpdate_plink_LOCO(n_nomissing, q, Y_LOCO, X, W_new, sparsekin, tau_new, fixtau, tol_tau, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);	
					
					arma::fvec temp2 = fit_LOCO["alpha"];
					alpha_new_LOCO = temp2;
					arma::fvec temp3 = fit_LOCO["eta"];
					eta_LOCO = temp3;

					if(verbose){
						Rcout << "Fixed-effect coefficients LOCO: " << alpha_new << std::endl;
					}

					// for gaussian family
					mu_LOCO = eta_LOCO;
					mu_eta_LOCO.fill(1.0);
					Y_LOCO = eta_LOCO + (y - mu_LOCO)/mu_eta_LOCO;

					// // check convergence
					arma::fvec check_vec_alpha = arma::abs(alpha_new_LOCO - alpha0_LOCO)/(arma::abs(alpha_new_LOCO) + arma::abs(alpha0_LOCO) + tol_coef);
					float check_max = check_vec_alpha.max();
					if( 2*check_max < tol_coef ){
						break;
					}
					
				} // end i < maxiter

				res_LOCO = y - mu_LOCO; 
				results_LOCO[chr-1] = List::create(Named("coefficients")=alpha_new_LOCO, Named("linear_predictors")=eta_LOCO, Named("fitted_values") = mu_LOCO, Named("Y")=Y_LOCO, Named("residuals")=res_LOCO, Named("cov")=fit_LOCO["cov"]);

			} // end if not na


		} // end for chr
		return List::create(Named("theta") = tau_new, Named("coefficients")=alpha_new, Named("linear_predictors")=eta, Named("fitted_values") = mu, Named("Y")=Y_new, Named("X")=X, Named("residuals")=res, Named("cov")=fit["cov"], Named("XSiX_inv_SiXt")=fit["XSiX_inv_SiXt"], Named("converged")=converged, Named("method") ="glmm_ai_multitau", Named("Sigma_iX")=Sigma_iX, Named("LOCO")=loco, Named("result_LOCO") = results_LOCO );

	} // end if loco

	
}



//' Variance component estimation
//' Uses Plink file for GRM construction and no sparse kinship matrix
//'
//' @title Iteratively estimate variance components tau in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param Y			column vector, working y
//' @param y			column vector, original y
//' @param X   			data matrix
//' @param W 			numeric vector.
//' @param alpha 		numeric vector.
//' @param maxiter 		integer, maximum iteration.
//' @param tol         	Numeric. Tolerance for estimation of variance components.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm. 
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm. 
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @param verbose		boolean. 
//' @param tol_tau		Numeric. Threshold for minimum variance components value. Default is 1e-5.
//' @param tol_coef    Numeric. Tolerance for convergence of coefficients. Default is 0.1
//' @return A List.
//' * theta: variance components estimates
//' * coefficients
//' * linear_predictors
//' * fitted_values
//' * Y
//' * X
//' * residuals
//' * cov
//' * XSiX_inv_SiXt
//' * converged
//' * method
// [[Rcpp::export]]
List getTau_plink_nok(int n_nomissing, arma::fvec tau, arma::fvec fixtau, arma::fvec& Y, arma::fvec& y, arma::fmat& X, arma::fvec alpha, arma::fvec W, int maxiter, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace, bool verbose, float tol_tau=1e-5, float tol_coef=0.1){

	int n =n_nomissing;
	// Kinships kinships(K, kmatfile);
	Function warning("warning");

	int q = 1; 
	
	arma::fvec W_new = W;
	// for iterative update
	arma::fvec alpha0, tau0;
	arma::fvec alpha_new = alpha, tau_new = tau;
	List fit(5);
	// arma::fmat cov(n,n);
	arma::fvec eta(n), mu(n), mu_eta(n);
	arma::fvec Y_new = Y; 

	int q2,i;
	arma::uvec idxtau = arma::find(fixtau<1);

	// ======== iterative update tau ===========
	for( i=0; i< maxiter; i++){
		if(verbose){
			Rcout << "------ Variance component estimation iteration " << i << " ------" << std::endl;
		}
		alpha0 = alpha_new;
		tau0 = tau_new;
		q2 = arma::sum(fixtau<1);

		// returns list with element "tau","alpha","eta","cov" (updated)
		fit = glmmaiUpdate_plink_nok(n_nomissing, q, Y, X, W_new, tau_new, fixtau, tol_tau, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);
		
		if(q2>0){
			arma::fvec temp1 = fit["tau"];
			tau_new = temp1;
			for(int j = 0; j<q+1; j++){
				if(fixtau(j)<1 && tau_new(j)< 1.01 * tol_tau){
					fixtau(j) = 1.0;
				}
			}		
		} // end if q2

		arma::fvec temp2 = fit["alpha"];
		alpha_new = temp2;
		arma::fvec temp3 = fit["eta"];
		eta = temp3;


		if(verbose){
			Rcout << "Variance component estimates: " << tau_new << std::endl;

		}

		// for gaussian family
		mu = eta;
		mu_eta.fill(1.0);
		Y_new = eta + (y - mu)/mu_eta;
		// W is always ones for gaussian family, no need to update

		// // check convergence
		arma::fvec check_vec_alpha = arma::abs(alpha_new - alpha0)/(arma::abs(alpha_new) + arma::abs(alpha0) + tol_coef);
		arma::fvec check_vec_tau = arma::abs(tau_new - tau0)/(arma::abs(tau_new) + arma::abs(tau0) + tol);
		
		if( (2*check_vec_tau.max() < tol) || ( 2*check_vec_alpha.max() < tol_coef) ){
			break;
		}

		float check_max = tau_new.max();
		if( check_max > 1.0/tol/tol) {
			warning("Large variance estimate observed in the iterations, model not converged...");
			i = maxiter;
			break;
		}
	} // end of iteration maxiter
	
	bool converged = (i < maxiter);
	arma::fvec res = y - mu; 

	arma::fmat Sigma_iX = pcgsolveMat_plink_nok(W, tau_new, X, "Jacobi", tol_pcg, maxiter_pcg);

	return List::create(Named("theta") = tau_new, Named("coefficients")=alpha_new, Named("linear_predictors")=eta, Named("fitted_values") = mu, Named("Y")=Y_new, Named("X")=X, Named("residuals")=res, Named("cov")=fit["cov"], Named("XSiX_inv_SiXt")=fit["XSiX_inv_SiXt"], Named("converged")=converged, Named("method") ="glmm_ai_multitau",Named("Sigma_iX")=Sigma_iX);
}