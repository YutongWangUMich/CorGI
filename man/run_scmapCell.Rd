% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmarking.R
\name{run_scmapCell}
\alias{run_scmapCell}
\title{Wrapper for \code{scmap::scmapCell}}
\usage{
run_scmapCell(query, ref, gene_set, param, ...)
}
\arguments{
\item{query}{the query dataset (SingleCellExperiment object)}

\item{ref}{the reference dataset (SingleCellExperiment object)}

\item{gene_set}{a list of characters representing the genes}

\item{param}{the number of nearest neighbors in scmapCell, an integer}

\item{...}{additional parameters passed on to \code{caret::confusionMatrix}. For example \code{dnn = c("Prediction","Truth")}.}
}
\value{
confusion matrix of the true labels in `ref` and the `scmapCluster` predicted labels
}
\description{
Wrapper for \code{scmap::scmapCell}
}
