% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pcgsolve_plink_LOCO}
\alias{pcgsolve_plink_LOCO}
\title{Solve for x in Sigma x = b using preconditioned conjugate gradient method.}
\usage{
pcgsolve_plink_LOCO(sparsekin, W, tau, b, preconditioner, tol, maxIter)
}
\arguments{
\item{sparsekin}{matrix of kinship coefficients (sparse).}

\item{W}{numeric vector.}

\item{tau}{numeric vector, initial values for variance components.}

\item{b}{vector, with same dimension as number of rows of A.}

\item{preconditioner}{string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.}

\item{tol}{numeric, threshold for convergence, default is \code{1e-6}.}

\item{maxIter}{numeric, maximum iteration, default is \code{1000}.}
}
\value{
A vector representing solution x.
}
\description{
Preconditioned conjugate gradient method LOCO (from cPCG)
}
\details{
Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm
}
