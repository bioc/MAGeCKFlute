% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IncorporateDepmap.R
\docType{methods}
\name{IncorporateDepmap}
\alias{IncorporateDepmap}
\title{Incorporate Depmap screen into analysis}
\usage{
IncorporateDepmap(dd, symbol = "id", cell_lines = NA,
  lineages = "All", na.rm = TRUE)
}
\arguments{
\item{dd}{A data frame.}

\item{symbol}{A character, specifying the column name of gene symbols in the data frame.}

\item{cell_lines}{A character vector, specifying the cell lines in Depmap to be considered.}

\item{lineages}{A character vector, specifying the lineages in Depmap to be considered.}

\item{na.rm}{Boolean, indicating whether removing NAs from the results.}
}
\value{
A data frame with Depmap column attached.
}
\description{
Incorporate Depmap screen into analysis
}
\examples{
dd.rra = ReadRRA(rra.gene_summary)
depmap_similarity = ResembleDepmap(dd.rra)
dd.rra = IncorporateDepmap(dd.rra, cell_lines=rownames(depmap_similarity)[1:3])
head(dd.rra)
}
\author{
Wubing Zhang
}