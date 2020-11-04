#ifndef PKG_GETTRACE_H
#define PKG_GETTRACE_H

arma::fvec getTrace(int n_nomissing, int q, arma::sp_mat& sparsekin, const std::vector<arma::fmat>& densekin, arma::fmat& Sigma, arma::fmat& Sigma_iX, arma::fmat& X, const arma::fvec& W_vec, arma::fmat& XSiX_inv,  int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace);

#endif