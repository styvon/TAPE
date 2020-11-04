#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <list>
#include <iterator>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::max


using namespace Rcpp; 
using namespace arma; 
#include "pcgsolve.h"
#include "getTrace.h"

//' Get score within a variance component estimation iteration
//'
//' @title Get score value in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param densekin 	list of dense covariance matrices read from kmatfile.
//' @param W 			numeric vector.
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param is_AI        boolean. Whether AI trace estimation will be used.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A list with following elements:
//' * YPAPY
//' * Trace
//' * PY
//' * AI
//' * cov = XSiX_inv
//' * tau
//' * alpha
//' * dtau
//' * eta
//' * XSiX_inv_SiXt
// [[Rcpp::export]]
List getScore(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const std::vector<arma::fmat>& densekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){

	Function warning("warning");
	// ==== variable declaration ====
	int n = n_nomissing;
	int q2 = arma::sum(fixtau==0);
	arma::uvec idxtau = arma::find(fixtau==0);

	int ncol_x = X.n_cols;
	arma::fvec tau0;

	arma::fmat Sigma_iX(n,ncol_x), Sigma_iXt(ncol_x,n), XSiX_inv(n,ncol_x);
	arma::fvec Sigma_iY(n);
	arma::fvec XmatVecTemp(n);
	arma::fvec PY(n);
	arma::vec PY_temp(n); // for sp mat format conversion
	
	arma::fmat APYmat(n,q+1); // store col vecs of APY for each var component
	arma::fvec PAPY_1(n), PAPY(n);
	// arma::fvec YPAPY(q+1); // for storing results from q+1 components
	arma::fmat AI(q+1,q+1);
	AI.zeros();
	arma::fvec YPAPY(q+1);
	YPAPY.zeros();
	arma::fvec Trace(q+1);
	Trace.zeros();

	tau0 = tau;
	
	// ===== get Sigma ====== 
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
			temp_dense = arma::conv_to<fmat>::from(densekin[i]);
			tempvec = (double)tau0(i+2)/arma::mean(arma::diagvec(temp_dense));
			tau0(i+2) = arma::conv_to<float>::from(tempvec);
			
			Sigma = Sigma + tau0(i+2)*temp_dense;
		}
	}

	Sigma = arma::symmatu(Sigma);	// force symmetric


	// ==== get computing components, alpha, eta ====
	for(int i = 0; i < ncol_x; i++){
        XmatVecTemp = X.col(i);
        // Sigma_iX.col(i) = pcgsolve(Sigma, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
        Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
    }
    Sigma_iXt = Sigma_iX.t();
    // Sigma_iY = pcgsolve(Sigma, Y, "Jacobi", tol_pcg, maxiter_pcg);
    Sigma_iY = arma::solve(Sigma, Y, arma::solve_opts::likely_sympd	);

    arma::fmat XSiX = arma::symmatu( Sigma_iXt * X) ;
    try {
    	// XSiX_inv = arma::inv(  X.t() * Sigma_iX  ); //=cov
    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
    }
    catch(const std::exception& e){
    	std::cout << "getScore: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
    	XSiX_inv = arma::pinv(   XSiX  ); //=cov
    }
    
    arma::fmat XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
    arma::fvec alpha = XSiX_inv_SiXt * Y;
    PY = Sigma_iY - Sigma_iX * alpha;

    arma::fvec diagp = tau0(0)/W;
    arma::fvec eta = Y - diagp % PY;
    // arma::fvec eta = Y - PY / (tau(0) / W);

    if(q2>0){
		// ==== get AI: 1/2====
		// for(int i=0; i<q2; i++){
	    for(int i=0; i<q+1; i++){
	    	// if(idxtau[i]==0){ // identity mat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	        	if(i==1){ // sparsekin
	        		// conversion for ops with sp_mat
	        		PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
	        		APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	        	}
	        	else { //densekin
	        		APYmat.col(i)= densekin[i-2]*PY;
	        	} // end if i==1 	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve(Sigma, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			// PAPY_1 = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
			// PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * PAPY_1));
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);
				if(j !=i){
					AI(j,i) = AI(i,j);
				}		
			}
			
			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i


		//====Calculate trace: 1/2====
		// use Hutchinson’s randomized trace estimator
		Trace = getTrace(n_nomissing, q, sparsekin, densekin, Sigma, Sigma_iX, X, W, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

		//==== update tau: 1/2 ====
		// tau0 = tau0 +  tau0 % tau0/float(n) % (YPAPY-Trace);
		for(int i=0; i<q2; i++){
			tau0( idxtau(i) ) = tau0( idxtau(i) ) +  tau0( idxtau(i) ) * tau0( idxtau(i) )/float(n) * (YPAPY( idxtau(i) )-Trace( idxtau(i) ));
		}// end for i in q2

		// ==== get AI: 2/2 ======
		for(int i = 0; i < ncol_x; i++){
	        XmatVecTemp = X.col(i);
	        // Sigma_iX.col(i) = pcgsolve(Sigma, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
	        Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
	    }
	    Sigma_iXt = Sigma_iX.t();
	    // Sigma_iY = pcgsolve(Sigma, Y, "Jacobi", tol_pcg, maxiter_pcg);
	    Sigma_iY = arma::solve(Sigma, Y, arma::solve_opts::likely_sympd	);
	    XSiX = arma::symmatu( Sigma_iXt * X) ;
	    try {
	    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
	    }
	    catch(const std::exception& e){
	    	std::cout << "getScore_plink: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
	    	XSiX_inv = arma::pinv(   XSiX  ); 
	    }   
	    XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
	    alpha = XSiX_inv_SiXt * Y;
	    PY = Sigma_iY - Sigma_iX * alpha;

	    for(int i=0; i<q+1; i++){
	    	// if(idxtau[i]==0){ // identity mat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	        	if(i==1){ // sparsekin
	        		// conversion for ops with sp_mat
	        		PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
	        		APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	        	}
	        	else { //densekin
	        		APYmat.col(i)= densekin[i-2]*PY;
	        	} // end if i==1 	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve(Sigma, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			// PAPY_1 = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
			// PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * PAPY_1));
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);
				if(j !=i){
					AI(j,i) = AI(i,j);
				}		
			}
			
			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i

	    //====Calculate trace: 2/2====
	    // use Hutchinson’s randomized trace estimator
	    Trace = getTrace(n_nomissing, q, sparsekin, densekin, Sigma, Sigma_iX, X, W, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

	    // ====== get dtau =======
		arma::fvec score_vec_update(q2);
		arma::fvec dtau(q+1), dtau_pre(q2);
		arma::fmat AI_update(q2,q2);

		for(int i=0; i<q2; i++){
			// fill AI_update
			for(int j=0; j<=i;j++){
				AI_update(i,j) = AI(idxtau(i),idxtau(j));
				if(j != i){
					AI_update(j,i) = AI_update(i,j);
				}
			} // end for j	
			// fill score_vec_update
			score_vec_update(i) = YPAPY(idxtau(i)) - Trace(idxtau(i));	
		} // end for i in q2
		try{
			// arma::fvec dtau = arma::solve(AI_mat, score_vec, arma::solve_opts::likely_sympd);
			dtau_pre = arma::solve(AI_update, score_vec_update, arma::solve_opts::allow_ugly);	
		}
		catch(std::runtime_error){
			std::cout << "getScore_plink: arma::solve(AI_mat, score_vec): AI seems singular, using less variant components matrix is suggested." << std::endl;
			dtau_pre.zeros();
		}

		// fill dtau using dtau_pre, padding 0
		int i2 = 0;
		for(int i=0; i<q+1; i++){
			if(fixtau(i)<1){ // not fixed
				dtau(i) = dtau_pre(i2);
				i2++;
			} else { // fixed
				dtau(i) = 0.0;
			}
		} // end for i

		//====Update tau: 2/2 in glmmaiUpdate ====
		// List out = List::create(Named("Sigma") = Sigma);
		List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=dtau, Named("eta") = eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);
		return(out);
	} // end if q2

	List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = NULL, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=NULL, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);

	return(out);

}

