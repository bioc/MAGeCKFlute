% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichedFilter.R
\docType{methods}
\name{EnrichedFilter}
\alias{EnrichedFilter}
\title{Simplify the enrichment results based on Jaccard index}
\usage{
EnrichedFilter(enrichment = enrichment, cutoff = 0.8)
}
\arguments{
\item{enrichment}{A data frame of enrichment result or an enrichResult object.}

\item{cutoff}{A numeric, specifying the cutoff of Jaccard index between two pathways.}
}
\value{
A data frame.
}
\description{
Simplify the enrichment results based on Jaccard index
}
\examples{
data(geneList, package = "DOSE")
\dontrun{
  enrichRes <- enrich.HGT(geneList, keytype = "entrez")
  EnrichedFilter(enrichRes)
}
}
\author{
Yihan Xiao
}
