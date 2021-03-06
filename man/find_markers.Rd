% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_markers.R
\name{find_markers}
\alias{find_markers}
\title{Remove genes (rows) with little information.
Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.}
\usage{
find_markers(X, genenames, cluster, tsne, top.n = 100, plot.n = 3)
}
\arguments{
\item{X}{log-transformed expression level matrix with cells at columns, genes at rows}

\item{genenames}{gene names}

\item{cluster}{: clustering result}

\item{tsne}{: tsne result}

\item{top.n}{: number of top genes to return for each cluster}

\item{plot.n}{: number of genes to plot for each cluster}
}
\value{
expression level matrix containing only highly variable genes
}
\description{
Remove genes (rows) with little information.
Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
}
