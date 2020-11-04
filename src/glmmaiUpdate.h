#ifndef PKG_GLMMAIUPDATE_H
#define PKG_GLMMAIUPDATE_H

List glmmaiUpdate(int n_nomissing, int q, arma::fvec Y, arma::fmat X, const arma::fvec& W, arma::sp_mat& sparsekin, const std::vector<arma::fmat>& densekin, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace);

#endif