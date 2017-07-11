

#' Extract the MAP pseudotime estimates from a \code{ouija_fit}
#'
#' @param oui An object of class \code{ouija_fit}.
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @export
#' @name map_pseudotime
#' 
#' @return MAP pseudotime vector of length N
#' 
#' @examples 
#' data(oui)
#' tmap <- map_pseudotime(oui)
map_pseudotime <- function(oui) {
  stopifnot(is(oui, "ouija_fit"))
  posterior.mode(mcmc(extract(oui$fit, "t")$t))
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
#' @importFrom coda HPDinterval
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
  hpd <- HPDinterval(mcmc(t_trace), prob = prob)
  return( hpd )
}


#' Predicted expression
#' 
#' The predicted (ie mean) expression from the fit.
#' 
#' @param oui An object of class \code{ouija_fit}.
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @return A matrix of the same dimension as \code{oui$Y} containing 
#' the predicted expression.
#' 
#' @export
#' @name rexprs
#' @examples 
#' data(oui)
#' pexp <- rexprs(oui)
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

sample_predicted_expression <- function(oui) {
  extr <- extract(oui$fit, pars = c("mu0", "k", "t0"))
  sample_ind <- sample(seq_len(nrow(extr$mu0)), 1)
  mu0 <- extr$mu0[sample_ind,]
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
  msg <- paste("A Bayesian non-linear factor analysis fit with\n",
               x$N, "cells and", x$G, "marker genes\n",
               "Inference type: ", itype)
  if(x$inference_type == "hmc") {
    msg <- paste(msg, "\nMCMC info:", x$iter, "iterations on", x$chains, "chains\n")
  }
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
#' @export
#' 
#' @examples 
#' tsigmoid(1, 1, 0.5, runif(10))
tsigmoid <- function(mu0, k, t0, t) {
  return( 2 * mu0 / (1 + exp(-k*(t - t0))))
}

transient <- function(mu0, p, b, t) {
  return( 2 * mu0 * exp(-10 * b * (t - p)^2) )
}

#' Create a confusion matrix
#' @export
confusion_matrix <- function(oui) {
  N <- oui$N
  pst_traces <- extract(oui$fit, "t")$t  
  confusion_mat <- matrix(NA, nrow = N, ncol = N)
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      confusion_mat[i,j] <- mean(pst_traces[,i] > pst_traces[,j])
      confusion_mat[j,i] <- 1 - confusion_mat[i,j]
    }
  }
  return(confusion_mat)
}

#' Order confusion matrix by pseudotime
#' @export
confusion_matrix_ordered <- function(oui, cmat = NULL) {
  if(is.null(cmat)) cmat <- confusion_matrix(oui)
  pst_order <- order(map_pseudotime(oui))
  cmat <- cmat[pst_order, pst_order]
  return(cmat)
}

#' Cluster the confusion matrix
#' @import mclust
#' @export
cluster_confusion <- function(cmat, n_clusters = 2:9) {
  diag(cmat) <- 0
  pca <- prcomp(cmat)
  pc1 <- pca$x[,1]
  mc <- Mclust(pc1, G = n_clusters)
  return(mc$classification)
}


