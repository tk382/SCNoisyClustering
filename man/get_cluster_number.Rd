% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cluster_number.R
\name{get_cluster_number}
\alias{get_cluster_number}
\title{Estimate the optimal number of clusters using eigengap}
\usage{
get_cluster_number(S)
}
\arguments{
\item{S}{estimated similarity matrix.}
}
\description{
Estimate the optimal number of clusters using eigengap
}
\examples{
#create positive definite symmetric matrix
X = matrix(rnorm(50), nrow = 10)
S= t(X) \%*\% X
getClustNum(S)
}
