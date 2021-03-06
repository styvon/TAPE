% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getTau_plink_nok}
\alias{getTau_plink_nok}
\title{Iteratively estimate variance components tau in GLMM}
\usage{
getTau_plink_nok(
  n_nomissing,
  tau,
  fixtau,
  Y,
  y,
  X,
  alpha,
  W,
  maxiter,
  tol,
  tol_pcg,
  maxiter_pcg,
  nrun_trace,
  cutoff_trace,
  verbose,
  tol_tau = 1e-05,
  tol_coef = 0.1
)
}
\arguments{
\item{n_nomissing}{integer, number of non-missing individuals as returned from setgeno()}

\item{tau}{numeric vector, initial values for variance components.}

\item{fixtau}{vector with elements 0/1, indicator for whether tau elements are fixed.}

\item{Y}{column vector, working y}

\item{y}{column vector, original y}

\item{X}{data matrix}

\item{alpha}{numeric vector.}

\item{W}{numeric vector.}

\item{maxiter}{integer, maximum iteration.}

\item{tol}{Numeric. Tolerance for estimation of variance components.}

\item{tol_pcg}{Numeric. Tolerance for the PCG algorithm.}

\item{maxiter_pcg}{Integer. Max number of iterations for the PCG algorithm.}

\item{nrun_trace}{Integer. Number of random vectors used for trace estimation.}

\item{cutoff_trace}{Numeric. Threshold for the coefficient of variation for trace estimation.}

\item{verbose}{boolean.}

\item{tol_tau}{Numeric. Threshold for minimum variance components value. Default is 1e-5.}

\item{tol_coef}{Numeric. Tolerance for convergence of coefficients. Default is 0.1}
}
\value{
A List.
* theta: variance components estimates
* coefficients
* linear_predictors
* fitted_values
* Y
* X
* residuals
* cov
* XSiX_inv_SiXt
* converged
* method
}
\description{
Variance component estimation
Uses Plink file for GRM construction and no sparse kinship matrix
}
