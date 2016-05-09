## Significance testing to weigh evidence that a set of genes
## are truly involved in pseudotemporal process
## Experimental!

#' Residual variance after fit
#' 
#' @param bm An object of class \code{bnlfa_fit}.
#' @param weights Cell specific weights.
#' 
#' @export
fit_rss <- function(bm, weights = NULL) {
  stopifnot(is(bm, "bnlfa_fit"))
  Z <- rexprs(bm)
  
  if(is.null(weights)) weights <- matrix(1, nrow = nrow(Z), ncol = ncol(Z))
  
  RSS <- colSums(weights * (bm$Y - Z)^2) # / (bm$N - 1)
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

#' F-test taking model complexity into account
reg_f_test <- function(bm) {
  #' From Wikipedia: Consider two models, 1 and 2, where model 1 is 'nested' 
  #' within model 2. Model 1 is the Restricted model, and Model 2 is the Unrestricted one. 
  #' That is, model 1 has p1 parameters, and model 2 has p2 parameters, 
  #' where p2 > p1, and for any choice of parameters in model 1, the same 
  #' regression curve can be achieved by some choice of the parameters of model 2. See
  #' https://en.wikipedia.org/wiki/F-test#Regression_problems
  #' 
  #' So for our model 1 is the null model (k=0) and 2 is the sigmoidal model, then
  #' 
  #' F = (n-p_2) * (RSS1 - RSS2) / {RSS2 * (p_2 - p_1)}
  #' 
  #' will have an F distribution with (p_2 - p_1, n - p_2) DOF.
  #' 
  #' Here p_2 = 4 and p_1 = 2
  weights <- 1 - prob_dropout(bm)
  
  RSS2 <- fit_rss(bm, weights)
  RSS1 <- null_rss(bm, weights)
  RSS_diff <- RSS1 - RSS2
  n <- dim(bm$Y)[1]
  p2 <- 4 
  p1 <- 2
  
  F_statistics <- (n - p2) * RSS_diff / ( RSS2 * (p2 - p1))
  lpvals <- pf(F_statistics, df1 = p2 - p1, df2 = n - p2, lower.tail = FALSE, log.p = TRUE)  
  p_val <- pchisq(-2 * sum(lpvals), df = 2 * length(lpvals), lower.tail = F)
  return(p_val)
}

#' @export
#' @importFrom rstan extract
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda mcmc
prob_dropout <- function(bm) {
  beta_map <- posterior.mode(mcmc(extract(bm$fit, "beta")$beta))

  pdr <- function(x, beta) 1 / exp(-(beta[1] + beta[2] * x))

  dropout_matrix <- apply(bm$Y, 2, pdr, beta_map) # gives probability of dropout
  return( dropout_matrix )
}

weighted_mean <- function(x, w) {
  return(sum(x * w) / sum(w))
}

null_rss <- function(bm, weights = NULL) {
  if(is.null(weights)) weights <- matrix(1, nrow = nrow(bm$Y), ncol = ncol(bm$Y))
  col_means <- sapply(1:ncol(bm$Y), function(i) weighted_mean(bm$Y[,i], weights[,i]))
  mean_matrix <- t(matrix(col_means, nrow = ncol(bm$Y), ncol = nrow(bm$Y)))
  RSS <- colSums(weights * (bm$Y - mean_matrix)^2)
  return(RSS)
}


