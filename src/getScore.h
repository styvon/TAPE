#ifndef PKG_GETSCORE_H
#define PKG_GETSCORE_H

List getScore(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const std::vector<arma::fmat>& densekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace);

#endif