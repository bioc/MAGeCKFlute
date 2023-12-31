% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichAB.R
\docType{methods}
\name{EnrichAB}
\alias{EnrichAB}
\title{Enrichment analysis for Positive and Negative selection genes}
\usage{
EnrichAB(
  data,
  enrich_method = "HGT",
  top = 10,
  limit = c(2, 100),
  filename = NULL,
  out.dir = ".",
  width = 6.5,
  height = 4,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{A data frame.}

\item{enrich_method}{One of "ORT" (Over-Representing Test) and "HGT" (HyperGemetric test).}

\item{top}{An integer, specifying the number of pathways to show.}

\item{limit}{A two-length vector, specifying the min and max size of pathways
for enrichent analysis.}

\item{filename}{Suffix of output file name.}

\item{out.dir}{Path to save plot to (combined with filename).}

\item{width}{As in ggsave.}

\item{height}{As in ggsave.}

\item{verbose}{Boolean}

\item{...}{Other available parameters in ggsave.}
}
\value{
A list containing enrichment results for each group genes. This list contains eight
items, which contain subitems of \code{gridPlot} and \code{enrichRes}.
}
\description{
Do enrichment analysis for selected genes, in which positive selection and negative selection
are termed as Positive and Negative
}
\author{
Wubing Zhang
}
