% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pcgsolveMat_plink_nok}
\alias{pcgsolveMat_plink_nok}
\title{Solve for x in Sigma x = B (matrix) using preconditioned conjugate gradient method.}
\usage{
pcgsolveMat_plink_nok(W, tau, bMat, preconditioner, tol, maxIter)
}
\arguments{
\item{W}{numeric vector.}

\item{tau}{numeric vector, initial values for variance components.}

\item{bMat}{vector, with same dimension as number of rows of A.}

\item{preconditioner}{string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.}

\item{tol}{numeric, threshold for convergence, default is \code{1e-6}.}

\item{maxIter}{numeric, maximum iteration, default is \code{1000}.}
}
\value{
A matrix representing solution Sigma_iB.
}
\description{
Preconditioned conjugate gradient method (from cPCG)
without sparse kinship matrix
}
\details{
Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm
}
