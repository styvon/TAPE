% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getDiagOfSigma}
\alias{getDiagOfSigma}
\title{Get diagonal elements from Sigma
Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm}
\usage{
getDiagOfSigma(sparsekin, W, tau)
}
\arguments{
\item{sparsekin}{matrix of kinship coefficients (sparse).}

\item{W}{numeric vector.}

\item{tau}{numeric vector, initial values for variance components.}
}
\value{
A vector with diagonal elements of Sigma
}
\description{
Get diagonal elements from Sigma
Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm
}
