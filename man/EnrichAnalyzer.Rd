% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichAnalyzer.R
\docType{methods}
\name{EnrichAnalyzer}
\alias{EnrichAnalyzer}
\alias{enrichment}
\title{Enrichment analysis}
\usage{
EnrichAnalyzer(
  geneList,
  keytype = "Symbol",
  type = "Pathway+GOBP",
  method = "HGT",
  organism = "hsa",
  pvalueCutoff = 1,
  limit = c(2, 100),
  universe = NULL,
  filter = FALSE,
  gmtpath = NULL,
  verbose = TRUE
)
}
\arguments{
\item{geneList}{A numeric vector with gene as names.}

\item{keytype}{"Entrez" or "Symbol".}

\item{type}{Molecular signatures for testing, available datasets include
Pathway (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF),
MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA), C2_CGP),
C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK)
and Complex (CORUM). Any combination of them are also accessible
(e.g. 'GOBP+GOMF+KEGG+REACTOME').}

\item{method}{One of "ORT"(Over-Representing Test), "GSEA"(Gene Set Enrichment Analysis), and "HGT"(HyperGemetric test).}

\item{organism}{'hsa' or 'mmu'.}

\item{pvalueCutoff}{FDR cutoff.}

\item{limit}{A two-length vector (default: c(2, 200)), specifying the minimal and
maximal size of gene sets for enrichent analysis.}

\item{universe}{A character vector, specifying the backgound genelist, default is whole genome.}

\item{filter}{Boolean, specifying whether filter out redundancies from the enrichment results.}

\item{gmtpath}{The path to customized gmt file.}

\item{verbose}{Boolean}
}
\value{
\code{enrichRes} is an enrichResult instance.
}
\description{
Enrichment analysis
}
\examples{
data(geneList, package = "DOSE")
\dontrun{
  keggA = EnrichAnalyzer(geneList[1:500], keytype = "entrez")
  head(keggA@result)
}
}
\seealso{
\code{\link{enrich.GSE}}

\code{\link{enrich.ORT}}

\code{\link{enrich.HGT}}

\code{\link[DOSE]{enrichResult-class}}
}
\author{
Wubing Zhang
}
