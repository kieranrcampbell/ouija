## Significance testing to weigh evidence that a set of genes
## are truly involved in pseudotemporal process
## Experimental!

#' @export
fit_rss <- function(bm) {
  stopifnot(is(bm, "bnlfa_fit"))
  Z <- rexprs(bm)
  RSS <- colSums((bm$Y - Z)^2) / (bm$N - 1)
  return(RSS)
}

#' @param ntrim Number of genes to retain (defaults to G / 2)
#' @export
#' @importFrom matrixStats colVars
fit_f_test <- function(bm, nkeep = NULL) {
  stopifnot(is(bm, "bnlfa_fit"))
  F_stat <- fit_rss(bm) / colVars(bm$Y)
  p_vals <- pf(F_stat, df1 = bm$N - 1, df2 = bm$N - 1)
  if(is.null(nkeep)) nkeep <- floor(0.5 * length(p_vals))
  trimmed_p_vals <- sort(p_vals, decreasing = TRUE)[seq_len(nkeep)]
  p_val <- pchisq(-2 * sum(log(trimmed_p_vals)), df = 2 * nkeep, lower.tail = F)
  return( p_val )
}

