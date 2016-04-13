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

#' @export
#' @importFrom matrixStats colVars
fit_f_test <- function(bm, prop_trim = 0.5) {
  stopifnot(is(bm, "bnlfa_fit"))
  F_stat <- fit_rss(bm) / colVars(bm$Y)
  p_vals <- pf(F_stat, df1 = bm$N - 1, df2 = bm$N - 1)
  np <- floor(length(p_vals) * prop_trim)
  trimmed_p_vals <- sort(p_vals, decreasing = TRUE)[seq_len(np)]
  p_val <- pchisq(-2 * sum(log(trimmed_p_vals)), df = 2 * np, lower.tail = F)
  return( p_val )
}

