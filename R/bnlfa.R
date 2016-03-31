#' Bayesian non-linear factor analysis for incorporating 
#' prior knowledge into single-cell trajectory learning
#' kieranc@well.ox.ac.uk


#' Fit a BNLFA object.
#' 
#' Fit a Bayesian non-linear factor analysis model given some single-cell
#' gene expression data.
#' 
#' This function takes either a \code{SCESet} object or an expression matrix
#' and returns a \code{bnlfa_fit} object including posterior traces for all
#' variables.
#' 
#' @param x Either an \code{SCESet} from \code{scater} or a
#' cell-by-gene (N by G) logged expression matrix
#' @param response The type of factor analysis, either \code{nonlinear} (default) 
#' or \code{linear} 
#' @param noise The pooling of the precision parameters, either \code{pool}, \code{partial-pool} 
#' (default) or \code{none}
#' @param prior The type of prior specification, either \code{normal} (default) or \code{sign}
#' @param model_mean_variance Logical. If \code{TRUE} then the variance as modelled as a function of the mean.
#' @param sign_bits A vector (length G) of sign bits
#' @param k_means G mean activation strength parameters
#' @param t0_means G mean activation time parameters
#' @param k_sd Optional standard deviations for k parameters
#' @param t0_sd Optional standard deviations for t0 parameters
#' @param ... Additional arguments to \code{rstan::sampling}
#' 
#' @import rstan
#' 
#' @export
#' 
#' @return An object of type \code{bnlfa_fit}
bnlfa <- function(x, response = c("nonlinear", "linear"),
                  noise = c("partial-pool", "pool", "none"),
                  prior = c("normal", "sign"),
                  model_mean_variance = TRUE,
                  sign_bits = NULL, k_means = NULL, t0_means = NULL,
                  k_sd = NULL, 
                  t0_sd = NULL,
                  lambda = 1e-2,
                  ...) {
  require(rstan) # for some reason this is required despite the @import rstan
  
  ## Find out what sort of model we're trying to fit
  response <- match.arg(response)
  noise <- match.arg(noise)
  prior <- match.arg(prior)
  
  if(!is.logical(model_mean_variance)) stop("Please specify mean_variance as logical")
  mean_variance <- as.integer(model_mean_variance)
  
  model_name <- paste(response, noise, prior, sep = "_")
  if(model_name != "nonlinear_partial-pool_normal" &&
     model_name != "nonlinear_none_normal") { # REMOVE WHEN MODELS IMPLEMENTED
    stop("Only nl-pp-n currently supported")
  }
  model_file <- paste0(model_name, ".stan")
  
  Y <- NULL
  if(is(x, "SCESet")) {
    ## convert to expression matrix Y  
    Y <- t(exprs(x))
  } else {
    Y <- x
  }
  if(!is(Y, "matrix")) {
    stop("x must either be an SCESet or matrix of gene expression values")
  }
  
  ## Now sanitize the input
  G <- ncol(Y) # number of genes
  N <- nrow(Y) # number of cells
  
  # we can fill in some values if they're null
  if(prior == "normal") {
    if(is.null(k_sd)) k_sd <- rep(1, G)
    if(is.null(t0_means)) t0_means <- rep(0.5, G)
    if(is.null(t0_sd)) t0_sd <- rep(1, G)
  }
  
  if(prior == "normal") {
    stopifnot(length(k_means) == G)
    stopifnot(length(k_sd) == G)
    if(response == "nonlinear") { # now we have t0 parameters
      stopifnot(length(t0_means) == G)
      stopifnot(length(t0_sd) == G)
    }
  } else {
    stopifnot(length(sign_bits) == G)
  }
  
  ## stan setup
  ## In a bit of glory for lazy programmers everywhere, rstan doesn't
  ## throw an error if we include stuff in data list that the model doesn't
  ## need. So we can fill this with junk for the other models. These can be left as NULL
  data <- list(Y = t(Y), G = G, N = N,
               k_means = k_means, k_sd = k_sd,
               t0_means = t0_means, t0_sd = t0_sd,
               lambda = lambda,
               mean_variance = mean_variance,
               sign_bits)
  
  stanfile <- system.file(model_file, package = "bnlfa")
  model <- stan_model(stanfile)
  
  ## manipulate stan defaults
  stanargs <- list(...)
  if(!('iter' %in% names(stanargs))) stanargs$iter <- 1e4
  if(!('chains' %in% names(stanargs))) stanargs$chains <- 1
  if(!('thin' %in% names(stanargs))) stanargs$thin <- 5
  stanargs$object <- model
  stanargs$data <- data
  
  ## call sampling
  fit <- do.call(sampling, stanargs)
  

  bm <- structure(list(fit = fit, G = G, N = N, Y = Y,
                       iter = stanargs$iter, chains = stanargs$chains,
                       thin = stanargs$thin, model_name = model_name), 
                  class = "bnlfa_fit")
  return(bm)
}

map_pseudotime <- function(bm) UseMethod("map_pseudotime")

#' Extract the MAP pseudotime values from a \code{bnlfa_fit}
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @export
#' 
#' @return MAP pseudotime vector of length N
map_pseudotime.bnlfa_fit <- function(bm) {
  posterior.mode(mcmc(extract(bm$fit, "t")$t))
}

#' Print a \code{bnlfa_fit}
#' 
#' @export
print.bnlfa_fit <- function(bm) {
  cat(paste("A Bayesian non-linear factor analysis fit of type", bm$model_name, "with\n"),
          paste(bm$N, "cells and", bm$G, "marker genes\n"),
          paste("MCMC info:", bm$iter, "iterations on", bm$chains, "chains"))
}

#' Plot a \code{bnlfa_fit}
#' 
#' @export
plot.bnlfa_fit <- function(bm, what = c("trace", "map", "diagnostic"), ...) {
  what <- match.arg(what)
  plt <- switch(what,
                trace = plot_bnlfa_fit_trace(bm, ...),
                map = plot_bnlfa_fit_map(bm, ...),
                diagnostic = plot_bnlfa_fit_diagnostics(bm, ...))
  return(plt)
}

#' Plot MCMC diagnostics.
#' 
#' Plot MCMC diagnostics (traceplot and autocorrelation) for a \code{bnlfa_fit} object
#' 
#' @param bm A \code{bnlfa_fit} object
#' @param nrow Number of rows. If 1, plots are side-by-side; if 2, plots are vertically aligned.
#' @export
#' @importFrom cowplot plot_grid
#' 
#' @return A \code{ggplot2} object
#' 
plot_bnlfa_fit_diagnostics <- function(bm, arrange = c("vertical", "horizontal")) {
  arrange <- match.arg(arrange)
  nrow <- switch(arrange,
                 vertical = 2,
                 horizontal = 1)
  plt <- cowplot::plot_grid(stan_trace(bm$fit, "lp__"), stan_ac(bm$fit, "lp__"), nrow = nrow)
  return(plt)
}

#' Plot heatmaps of gene expression changes
#' 
#' @importFrom rstan extract
#' @importFrom cowplot plot_grid
#' @importFrom viridis scale_fill_viridis
#' @importFrom reshape2 melt
#' 
#' 
#' @param samples Number of posterior pseudotime samples to use
#' @export
#' 
plot_bnlfa_fit_trace <- function(bm, samples = 50, genes = 1:min(bm$G, 5),
                                 output = c("grid", "plotlist"), 
                                 show_legend = FALSE, ...) {
  output <- match.arg(output)
  
  ttrace <- extract(bm$fit, "t")$t
  to_sample <- sample(seq_len(min(nrow(ttrace), samples)))
  ttrace <- ttrace[to_sample, ]
  cell_orders <- apply(ttrace, 1, order)
  # apply over genes
  plts <- lapply(genes, function(g) {
    yg <- bm$Y[,g]
    Xg <- apply(cell_orders, 2, function(or) yg[or])
    Xg <- data.frame(Xg)
    names(Xg) <- 1:ncol(Xg)
    Xg$x <- 1:nrow(Xg)
    Xm <- melt(Xg, id.vars = "x", variable.name = "y", value.name = "Expression")
    
    plt <- ggplot(Xm, aes(x = x, y = y, fill = Expression)) + geom_tile() +
      xlab("Cell") + ylab(expression("Pseudotime\n sample")) + 
      scale_fill_viridis() + 
      theme(axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank())
    if(!show_legend) plt <- plt + theme(legend.position = "none")
    return( plt )
  })

  if(output == "grid") {
    return(plot_grid(plotlist = plts, ...))
  } else {
    return( plts )
  }
}

#' Plot gene expression as a function of MAP pseudotime
#' 
#' @export
plot_bnlfa_fit_map <- function(bm) {
  tmap <- map_pseudotime(bm)
  Y <- bm$Y
  dy <- data.frame(Y, pseudotime = tmap)
  dm <- reshape2::melt(dy, id.vars = "pseudotime", 
                               variable.name = "gene", 
                               value.name = "expression")
  plt <- ggplot(dm, aes(x = pseudotime, y = expression)) + geom_point() +
    stat_smooth(color = 'red') + facet_wrap(~ gene, scales = "free_y") +
    xlab("MAP pseudotime")
  return(plt)
}

#' Synthetic gene expression matrix
#' 
#' A matrix containing some synthetic gene expression data for 100 cells and 3 genes
#' 
"synth_gex"

#' Synthetic gene pseudotimes
#' 
#' A vector with the 'true' pseudotimes for the synthetic gene expression data in \code{synth_gex}
"true_pst"