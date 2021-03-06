% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getTrace_plink}
\alias{getTrace_plink}
\title{Get trace using Hutchinson’s randomized trace estimator}
\usage{
getTrace_plink(
  n_nomissing,
  q,
  sparsekin,
  W,
  tau,
  Sigma_iX,
  XSiX_inv,
  nrun_trace,
  maxiter_pcg,
  tol_pcg,
  cutoff_trace
)
}
\arguments{
\item{n_nomissing}{integer}

\item{q}{integer, total number of sparse and dense correlation matrices}

\item{sparsekin}{matrix of kinship coefficients (sparse).}

\item{W}{numeric vector.}

\item{tau}{numeric vector, initial values for variance components.}

\item{Sigma_iX}{matrix}

\item{XSiX_inv}{matrix, covariance}

\item{nrun_trace}{Integer. Number of random vectors used for trace estimation.}

\item{maxiter_pcg}{Integer. Max number of iterations for the PCG algorithm.}

\item{tol_pcg}{Numeric. Tolerance for the PCG algorithm.}

\item{cutoff_trace}{Numeric. Threshold for the coefficient of variation for trace estimation.}
}
\description{
Get trace using Hutchinson’s randomized trace estimator
}
