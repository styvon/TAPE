#ifndef PKG_PCGUTIL_H
#define PKG_PCGUTIL_H

arma::fmat icc(arma::fmat A);
void set_seed(unsigned int seed);
Rcpp::NumericVector nb(int n);
float calCV(arma::fvec& xVec);

#endif