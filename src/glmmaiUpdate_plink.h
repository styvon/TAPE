#ifndef PKG_GLMMAIUPDATE_PLINK_H
#define PKG_GLMMAIUPDATE_PLINK_H

List glmmaiUpdate_plink(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::sp_mat& sparsekin, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace);
List glmmaiUpdate_plink_nok(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace);
List glmmaiUpdate_plink_LOCO(int n_nomissing, int q, arma::fvec& Y, arma::fmat X, const arma::fvec& W, arma::sp_mat& sparsekin, arma::fvec tau, arma::fvec fixtau, float tol, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace);


#endif