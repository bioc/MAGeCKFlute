#' Simplify the enrichment results based on Jaccard index
#'
#' @docType methods
#' @name EnrichedFilter
#' @rdname EnrichedFilter
#'
#' @param enrichment A data frame of enrichment result or an enrichResult object.
#' @param cutoff A numeric, specifying the cutoff of Jaccard index between two pathways.
#'
#' @return A data frame.
#'
#' @author Yihan Xiao
#' @examples
#' data(geneList, package = "DOSE")
#' \dontrun{
#'   enrichRes <- enrich.HGT(geneList, keytype = "entrez")
#'   EnrichedFilter(enrichRes)
#' }
#' @export

EnrichedFilter <- function(enrichment = enrichment, cutoff = 0.8){
  if(is(enrichment, "enrichResult")) enrichment = enrichment@result
  if(is(enrichment, "gseaResult")) enrichment = enrichment@result

  if(nrow(enrichment)<3) return(enrichment)
  enrichment = enrichment[order(enrichment$pvalue, -abs(enrichment$NES)), ]
  genelist = strsplit(enrichment$geneID, "/")
  names(genelist) = enrichment$ID
  # Jaccard Index
  tmp1 = crossprod(table(stack(genelist)))
  tmp2 = outer(lengths(genelist), lengths(genelist), "+")
  ijc = tmp1 / (tmp2 - tmp1)
  diag(ijc) = 0
  idx <- which(ijc>cutoff, arr.ind = TRUE)
  colnames(idx) = c("row", "col")
  idx = unlist(apply(idx, 1, max))
  enrichment <- enrichment[setdiff(1:nrow(enrichment), idx), ]

  # The second round simplification
  genelist = strsplit(enrichment$geneID, "/")
  names(genelist) = enrichment$ID
  tmp1 = crossprod(table(stack(genelist)))
  ijc2 = tmp1 / lengths(genelist)
  diag(ijc2) = 0
  idx <- which(ijc2>cutoff, arr.ind = TRUE)
  colnames(idx) = c("row", "col")
  idx = unlist(apply(idx, 1, max))
  enrichment <- enrichment[setdiff(1:nrow(enrichment), idx), ]
  return(enrichment)
}
