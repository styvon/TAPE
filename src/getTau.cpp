#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
#include <list>
#include <iterator>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace Rcpp; 
#include "glmmaiUpdate.h"
#include "readMat.h"

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


//' Variance component estimation
//'
//' @title Iteratively estimate variance components tau in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param K 			matrix of kinship coefficients (sparse).
//' @param kmatfile 	string vector, path to covariance matrices.
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
//' @param tol_tau		Numeric. Threshold for minimum variance components. Default is 1e-5.
//' @return A List.
//' * theta: variance components estimates
//' * coefficients
//' * linear_predictors
//' * fitted_values
//' * Y
//' * X
//' * residuals
//' * cov
//' * converged
//' * method
// [[Rcpp::export]]
List getTau(int n_nomissing, arma::sp_mat& sparsekin, StringVector kmatfile, arma::fvec tau, arma::fvec fixtau, arma::fvec& Y, arma::fvec& y, arma::fmat& X, arma::fvec alpha, arma::fvec W, int maxiter, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace, bool verbose, float tol_tau=1e-5){

	int n = n_nomissing;
	int ncol_x = X.n_cols;
	// Kinships kinships(K, kmatfile);
	Function warning("warning");
	// for store cov mats
	// arma::sp_fmat sparsekin(K);
	std::vector<arma::fmat> densekin;
	int q_de, q, i;
	
	arma::fvec W_new = W;
	// for iterative update
	arma::fvec alpha0, tau0;
	arma::fvec alpha_new = alpha, tau_new = tau;
	List fit(5);
	arma::fmat cov(n,n);
	arma::fvec eta(n), mu(n), mu_eta(n);
	arma::fvec Y_new = Y; 

	// new for dependence on fixtau
	int q2 = arma::sum(fixtau<1);
	arma::uvec idxtau = arma::find(fixtau<1);
	arma::fvec fixtau_old = fixtau;
	arma::fvec fixtau_new;

	arma::fmat temp;
	std::string tempfile;
	if(kmatfile(0)==""){
		// ======== store covariance matrices ===========
		q_de = 0;
		q = q_de +1;
		densekin.clear(); // set to empty

	}else{ // if has dense cov matrix in kmatfile
		// ======== store covariance matrices ===========
		q_de = kmatfile.length();
		q = q_de +1;
		
		for (i = 0; i< q_de; i++){
			tempfile = kmatfile(i);
			// .rel: "\t" delimited, else assume " " delimited
			if (hasEnding(tempfile,"rel")){
				temp = readMat(tempfile,"\t");
				temp.resize(n,n);
			}else{
				// temp.load(tempfile, arma::raw_ascii);
				// temp = readMat(tempfile,",");
				temp.load(tempfile, arma::csv_ascii);
				temp.resize(n,n);
			}	
			densekin.push_back(temp);
		}

	} // end if kmatfile=="" else
	// ======== iterative update tau ===========
	for( i=0; i< maxiter; i++){
		if(verbose){
			Rcout << "Variance component estimation iteration " << i << ":" << std::endl;
		}

		alpha0 = alpha_new;
		tau0 = tau_new;
		q2 = arma::sum(fixtau<1);

		// returns list with element "tau","alpha","eta","cov" (updated)
		fit = glmmaiUpdate(n_nomissing, q, Y_new, X, W_new, sparsekin, densekin, tau_new, fixtau, tol, tol_pcg, maxiter_pcg, nrun_trace, cutoff_trace);
		
		if(q2>0){
			arma::fvec temp1 = fit["tau"];
			tau_new = temp1;
			for(int j=0; j < q+1; j++){
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
			Rcout << "Fixed-effect coefficients: " << alpha_new << std::endl;

		}

		// for gaussian family
		mu = eta;
		mu_eta.fill(1.0);
		Y_new = eta + (y - mu)/mu_eta;
		// W is always ones for gaussian family, no need to update

		// check convergence
		arma::fvec check_vec_alpha = arma::abs(alpha_new - alpha0)/(arma::abs(alpha_new) + arma::abs(alpha0) + tol);
		arma::fvec check_vec_tau = arma::abs(tau_new - tau0)/(arma::abs(tau_new) + arma::abs(tau0) + tol);
		// arma::fvec check_vec_alpha = arma::abs(alpha_new - alpha0)/(arma::abs(alpha_new) + arma::abs(alpha0) );
		// arma::fvec check_vec_tau = arma::abs(tau_new - tau0)/(arma::abs(tau_new) + arma::abs(tau0) );

		float check_max = std::max(check_vec_alpha.max(), check_vec_tau.max());
		if( 2* check_max < tol ){
			break;
		}

		check_max = tau_new.max();
		if( check_max > 1.0/tol/tol) {
			warning("Large variance estimate observed in the iterations, model not converged...");
			i = maxiter;
			break;
		}
	} // end of iteration maxiter
	
	bool converged = (i < maxiter);
	arma::fvec res = y - mu; 

	// ===== get Sigma, Sigma_iX ====== 
	arma::fmat Wdiag = arma::diagmat(1/W); // Wdiag is 1/W in paper

	// convert formula with sp_mat to mat (no direct method provided by arma)
	arma::vec tempvec;
	tempvec = (double)tau0(1)/arma::mean(arma::diagvec(sparsekin));
	tau0(1) = arma::conv_to<float>::from(tempvec);
	arma::mat tempmat = 0.0 + (double)tau0(1) * sparsekin; 	
	arma::fmat Sigma = tau0(0) * Wdiag + arma::conv_to<arma::fmat>::from(tempmat);

	arma::fmat temp_dense(n,n);
	if(!densekin.empty()){ // if densekin not empty
		for(int i = 0; i< densekin.size(); i++){
			temp_dense = arma::conv_to<arma::fmat>::from(densekin[i]);
			tempvec = (double)tau0(i+2)/arma::mean(arma::diagvec(temp_dense));
			tau0(i+2) = arma::conv_to<float>::from(tempvec);		
			Sigma = Sigma + tau0(i+2)*temp_dense;
		}
	}
	Sigma = arma::symmatu(Sigma);	// force symmetric

	arma::fmat Sigma_iX(n,ncol_x);
	arma::fvec XmatVecTemp(n);
	for(int i = 0; i < ncol_x; i++){
        XmatVecTemp = X.col(i);
        // Sigma_iX.col(i) = pcgsolve(Sigma, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
        Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
    }

	return List::create(Named("theta") = tau_new, Named("coefficients")=alpha_new, Named("linear_predictors")=eta, Named("fitted_values") = mu, Named("Y")=Y, Named("X")=X, Named("residuals")=res, Named("cov")=fit["cov"], Named("XSiX_inv_SiXt")=fit["XSiX_inv_SiXt"], Named("converged")=converged, Named("method") ="glmm_ai_multitau", Named("Sigma_iX")=Sigma_iX);
}