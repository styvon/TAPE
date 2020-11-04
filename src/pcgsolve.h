#ifndef PKG_PCGSOLVE_H
#define PKG_PCGSOLVE_H

#include "pcgUtil.h"
arma::fvec pcgsolve(arma::fmat & A, arma::fvec & b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);

#endif