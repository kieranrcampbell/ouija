

#' Extract the MAP estimates from a \code{ouija_fit}
#' 
#' Extract the MAP (maximum a posteriori, aka posterior mode) 
#' estimates of the pseudotimes, switch times, peak times
#' or switch strengths from a \code{ouija_fit}.
#'
#' @param oui An object of class \code{ouija_fit}.
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @export
#' @name map
#' @importFrom methods is
#' 
#' @return MAP estimates of the specified quantity.
#' 
#' @examples 
#' data(oui)
#' tmap <- map_pseudotime(oui)
#' t0map <- switch_times(oui)
#' pmap <- peak_times(oui)
#' kmap <- switch_strengths(oui)
map_pseudotime <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  posterior.mode(mcmc(extract(oui$fit, "t")$t))
}

#' @name map
#' @export
switch_times <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  if(!any(oui$response_type == "switch")) return(NULL)
  posterior.mode(mcmc(extract(oui$fit, "t0")$t0))
}

#' @name map
#' @export
peak_times <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  if(!any(oui$response_type == "transient")) return(NULL)
  posterior.mode(mcmc(extract(oui$fit, "p")$p))
}

#' @name map
#' @export
switch_strengths <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  if(!any(oui$response_type == "switch")) return(NULL)
  posterior.mode(mcmc(extract(oui$fit, "k")$k))
}


#' Pseudotime errors
#' 
#' Returns the highest probability credible intervals for each cell.
#' @param oui An object of class \code{ouija_fit}
#' @param prob The probability for the credible interval. Default is 0.95.
#' 
#' @name pseudotime_error
#' 
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' 
#' @return A matrix with two columns, where the first column gives the lower
#' interval and the second gives the upper interval for the pseudotimes of
#' each cell.
#' 
#' @examples 
#' data(oui)
#' pst_err <- pseudotime_error(oui)
pseudotime_error <- function(oui, prob = 0.95) {
  t_trace <- extract(oui$fit, "t")$t
  hpd <- coda::HPDinterval(mcmc(t_trace), prob = prob)
  return( hpd )
}


#' Predicted expression (mean and sampling)
#' 
#' Extract the posterior mean of the predicted expression (ie the sigmoid
#' or transient function), or a sample from the posterior.
#' 
#' @param oui An object of class \code{ouija_fit}.
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @return A matrix of the same dimension as \code{oui$Y} containing 
#' the MAP mean expression or a posterior sample of the mean expression.
#' 
#' @export
#' @name predexprs
#' @examples 
#' data(oui)
#' pexp <- predicted_expression(oui)
#' sample_pexp <- sample_predicted_expression(oui)
predicted_expression <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  # Switch parameters
  extr_switch <- extract(oui$fit, pars = c("mu0_switch", "k", "t0"))
  mu0_map <- posterior.mode(mcmc(extr_switch$mu0_switch))
  k_map <- posterior.mode(mcmc(extr_switch$k))
  t0_map <- posterior.mode(mcmc(extr_switch$t0))
  sig_map <- mapply(tsigmoid, mu0_map, k_map, t0_map, 
                   MoreArgs = list(t = map_pseudotime(oui)))
  
  # Transient parameters
  if(any(oui$response_type == "transient")) {
    extr_transient <- extract(oui$fit, c("mu0_transient", "p", "b"))
    mu0_transient_map <- posterior.mode(mcmc(extr_transient$mu0_transient))
    p_map <- posterior.mode(mcmc(extr_transient$p))
    b_map <- posterior.mode(mcmc(extr_transient$b))
    t_map <- mapply(transient, mu0_transient_map, p_map, b_map,
                    MoreArgs = list(t = map_pseudotime(oui)))
  }
  
  pred_expr <- matrix(NA, nrow = nrow(oui$Y), ncol = ncol(oui$Y))
  colnames(pred_expr) <- colnames(oui$Y)
  pred_expr[, oui$response_type == "switch"] <- sig_map
  if(any(oui$response_type == "transient")) {
    pred_expr[, oui$response_type == "transient"] <- t_map
  }
  return(pred_expr)
}

#' @name predexprs
#' @export
sample_predicted_expression <- function(oui) {
  extr <- extract(oui$fit, pars = c("mu0_switch", "k", "t0"))
  sample_ind <- sample(seq_len(nrow(extr$mu0_switch)), 1)
  mu0 <- extr$mu0_switch[sample_ind,]
  k <- extr$k[sample_ind,]
  t0 <- extr$t0[sample_ind,]
  sig_map <- mapply(tsigmoid, mu0, k, t0, 
                    MoreArgs = list(t = map_pseudotime(oui)))
  colnames(sig_map) <- colnames(oui$Y)
  return(sig_map)
}


#' Print a \code{ouija_fit}
#' 
#' @param x An object of class \code{ouija_fit}.
#' @param ... Additional arguments.
#' 
#' @method print ouija_fit
#' @export
#' @return A character string representation of \code{x}
#' @examples 
#' data(oui)
#' print(oui)
print.ouija_fit <- function(x, ...) {
  itype <- switch(x$inference_type,
                  hmc = "Hamiltonian Monte Carlo",
                  vb = "Variational Bayes")
  msg <- paste("A Ouija fit with",
               x$N, "cells and", x$G, "marker genes",
               "\nInference type: ", itype)
  if(x$inference_type == "hmc") {
    msg <- paste(msg, "\nMCMC info:", x$iter, "iterations on", x$chains, "chains")
  }
  n_switch <- sum(x$response_type == "switch")
  n_transient <- sum(x$response_type == "transient")
  msg <- paste(msg, "\n(Gene behaviour) Switch/transient:",n_switch, "/", n_transient)
  cat(msg)
}


#' Generic function to return sigmoid whenever needed
#' 
#' @param mu0 The 'average' expression parameter
#' @param k The 'activation strength' parameter
#' @param t0 The 'activation time' parameter
#' @param t The latent pseudotimes
#' 
#' @return A numeric vector of length \code{length(t)} that
#' is the sigmoid function applied to \code{t} using the parameters
#' \code{mu0}, \code{k} and \code{t0}.
#' 
#' @keywords internal
#'
tsigmoid <- function(mu0, k, t0, t) {
  return( 2 * mu0 / (1 + exp(-k*(t - t0))))
}


#' Generic function to model transient expression whenever needed
#' 
#' @param mu0 Half peak expression
#' @param p Peak time
#' @param b Peak length
#' @param t Pseudotime
#' 
#' @keywords internal
#' @return The transient (rbf) function evaluated at input coordinates
transient <- function(mu0, p, b, t) {
  return( 2 * mu0 * exp(-10 * b * (t - p)^2) )
}

#' Create, order, and plot a consistency matrix
#' 
#' The consistency matrix is the N-by-N (for N cells) matrix where the entry in the
#' i-th row and j-th column describes the empirical probability that cell i is ordered
#' after cell j. These functions allow the calculation, ordering (by posterior pseudotime)
#' and plotting of the consistency matrix. See \code{cluster_consistency} to cluster
#' the cells using Gaussian Mixture Modelling.
#' 
#' The consistency matrix is useful for clustering the pseudotime trajectory into discrete
#' stages, as if there are regions where the empirical ordering probability is around 0.5
#' then there is large uncertainty as to the ordering, implying the cells are of roughly
#' the same type. If the probability is closer to 0 or 1 then there is less uncertainty
#' and the cells are probably undergoing a smooth transition.
#' 
#' @param oui A \code{ouija_fit} object
#' @seealso cluster_consistency
#' 
#' @name consistency
#' @export
#' @examples 
#' cmat <- consistency_matrix(oui)
#' 
#' @return A cell-by-cell consistency matrix (ordered by posterior pseudotime
#' if \code{consistency_matrix_ordered}), or a \code{ggplot2} plot object
#' displaying the ordered consistency matrix (if \code{plot_consistency_matrix})
consistency_matrix <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  N <- oui$N
  pst_traces <- extract(oui$fit, "t")$t  
  consistency_mat <- matrix(NA, nrow = N, ncol = N)
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      consistency_mat[i,j] <- mean(pst_traces[,i] > pst_traces[,j])
      consistency_mat[j,i] <- 1 - consistency_mat[i,j]
    }
  }
  return(consistency_mat)
}

#' @param cmat A pre-computed consistency matrix returned by \code{consistency_matrix(oui)}
#' @name consistency
#' @export
#' @examples
#' cmo <- consistency_matrix_ordered(oui, cmat)
consistency_matrix_ordered <- function(oui, cmat = NULL) {
  if(is.null(cmat)) cmat <- consistency_matrix(oui)
  pst_order <- order(map_pseudotime(oui))
  cmat <- cmat[pst_order, pst_order]
  return(cmat)
}

#' Cluster the consistency matrix
#' 
#' Clusters the consistency matrix, using Gaussian Mixture Modelling
#' on its first principal component. The number of clusters is chosen
#' such that the BIC is maximise
#' 
#' @param cmat A consistency matrix returned by \code{consistency_matrix}
#' @param n_clusters The number of clusters for which to calculate the BIC
#' 
#' @import mclust
#' @import dplyr
#' @export
#' @seealso consistency_matrix
#' @examples 
#' data(oui)
#' cmat <- consistency_matrix(oui)
#' clusters <- cluster_consistency(cmat)
#' 
#' @return A numeric vector indicating to which cluster each cell has been assigned
cluster_consistency <- function(cmat, n_clusters = 2:9) {
  diag(cmat) <- 0
  pca <- prcomp(cmat)
  pc1 <- pca$x[,1]
  mc <- Mclust(pc1, G = n_clusters)
  return(mc$classification)
}

#' Get regulation df
#' @importFrom rstan extract
#' @import dplyr
#' @keywords internal
#' @return A \code{data.frame} containing the posterior traces
#' for all pairwise differences in regulation timing
regulation_df <- function(oui) {
  y_index <- . <- NULL
  response_type <- oui$response_type
  G_switch <- sum(response_type == "switch")
  G_transient <- sum(response_type == "transient")
  gene_names <- colnames(oui$Y)
  
  switch_times <- extract(oui$fit, "t0")$t0
  peak_times <- extract(oui$fit, "p")$p
  
  switch_map <- data.frame(
    switch_index = seq_len(G_switch),
    y_index = which(response_type == "switch")
  )
  
  transient_map <- data.frame(
    transient_index = seq_len(G_transient),
    y_index = which(response_type == "transient")
  )
  
  dfs <- list()
  
  for(i in 1:(oui$G-1)) {
    for(j in (i+1):oui$G) {
      gene_i <- gene_names[i]
      rtype_i <- response_type[i]
      gene_j <- gene_names[j]
      rtype_j <- response_type[j]
      param_i <- param_j <- NULL
      
      if(rtype_i == "switch") param_i <- switch_times[,dplyr::filter(switch_map, y_index == i) %>% .$switch_index]
      if(rtype_j == "switch") param_j <- switch_times[,dplyr::filter(switch_map, y_index == j) %>% .$switch_index]
      if(rtype_i == "transient") param_i <- peak_times[,dplyr::filter(transient_map, y_index == i) %>% .$transient_index]
      if(rtype_j == "transient") param_j <- peak_times[,dplyr::filter(transient_map, y_index == j) %>% .$transient_index]
      
      df <- data_frame(gene_i = gene_i, rtype_i = rtype_i,
                       gene_j = gene_j, rtype_j = rtype_j,
                       param_diffs = param_i - param_j)
      dfs[[length(dfs) + 1]] <- df
    }
  }
  df_all <- bind_rows(dfs) %>% as_data_frame()
  df_all <- mutate(df_all,
                   label = paste(gene_i, "-", gene_j))
  return(df_all)
}


#' Differences in regulation timing
#' 
#' 
#' Computes a data frame with a row for every pairwise comparison
#' between genes, with a column for the posterior mean difference 
#' in regulation timings, along with the lower and upper 95% 
#' highest posterior density credible intervals, and a logical
#' indicating whether the differences in regulation timings
#' are significant.
#' 
#' @param oui An object of type \code{ouija_fit}
#' @name genereg
#' @export
#' @importFrom coda mcmc
#' @examples 
#' data(oui)
#' genereg <- gene_regulation(oui)
#' 
#' @return A \code{data.frame} as described above
gene_regulation <- function(oui) {
  label <- param_diffs <- gene_i <- gene_j <- gene_A <- gene_B <- NULL
  rtype_i <- rtype_j <- NULL
  
  reg_df <- regulation_df(oui)
  reg_df <- dplyr::select(reg_df, -rtype_i, -rtype_j)
  
  get_95 <- function(v, i) coda::HPDinterval(mcmc(v))[1,i]
  reg_df <- reg_df %>% 
    group_by(label, gene_i, gene_j) %>% 
    dplyr::summarise(mean_difference = mean(param_diffs),
              lower_95 = get_95(param_diffs, 1),
              upper_95 = get_95(param_diffs, 2))
  reg_df <- mutate(reg_df, signif = FALSE)
  reg_df$signif[reg_df$mean_difference < 0 & reg_df$upper_95 < 0] <- TRUE
  reg_df$signif[reg_df$mean_difference > 0 & reg_df$lower_95 > 0] <- TRUE
  reg_df <- dplyr::rename(reg_df, gene_A = gene_i, gene_B = gene_j, 
                          significant = signif)
  reg_df
}

