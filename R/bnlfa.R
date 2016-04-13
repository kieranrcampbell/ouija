## Bayesian non-linear factor analysis for incorporating 
## prior knowledge into single-cell trajectory learning
## kieranc@well.ox.ac.uk


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
#' cell-by-gene (N by G) matrix of non-negative values representing gene expression.
#' log2(TPM + 1) is recommended.
#' @param response The type of factor analysis, either \code{nonlinear} (default) 
#' or \code{linear} 
#' @param noise_pooling The pooling of the precision parameters, either \code{pool}, \code{partial-pool} 
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
                  noise_pooling = c("partial-pool", "pool", "none"),
                  prior = c("normal", "sign"),
                  init = "random",
                  model_mean_variance = TRUE,
                  sign_bits = NULL, k_means = NULL, t0_means = NULL,
                  k_sd = NULL, 
                  t0_sd = NULL,
                  lambda = 1,
                  ...) {
  require(rstan) # for some reason this is required despite the @import rstan
  
  ## Find out what sort of model we're trying to fit
  response <- match.arg(response)
  noise <- match.arg(noise_pooling)
  prior <- match.arg(prior)
  
  if(!is.logical(model_mean_variance)) stop("Please specify mean_variance as logical")
  mean_variance <- as.integer(model_mean_variance)
  
  model_name <- paste(response, noise, prior, sep = "_")
  if(model_name != "nonlinear_partial-pool_normal" &&
     model_name != "nonlinear_none_normal" &&
     model_name != "nonlinear_partial-pool_sign") { # REMOVE WHEN MODELS IMPLEMENTED
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
    if(is.null(k_means)) k_means = rep(0, G)
    if(is.null(k_sd)) k_sd <- rep(1, G)
    if(is.null(t0_means)) t0_means <- rep(0.5, G) ## change if constrained
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
               sign_bits = sign_bits)
  
  stanfile <- system.file(model_file, package = "bnlfa")
  model <- stan_model(stanfile)
  
  ## manipulate stan defaults
  stanargs <- list(...)
  if(!('iter' %in% names(stanargs))) stanargs$iter <- 1e4
  if(!('warmup' %in% names(stanargs))) stanargs$warmup <- stanargs$iter / 2
  if(!('chains' %in% names(stanargs))) stanargs$chains <- 1
  if(!('thin' %in% names(stanargs))) {
    # always specify thin so that approximately 1000 samples are returned
    stanargs$thin <- ceiling((stanargs$iter - stanargs$warmup) / 1000)
  }
  stanargs$object <- model
  stanargs$data <- data
  stanargs$init <- init
  
  ## call sampling
  fit <- do.call(sampling, stanargs)
  

  bm <- structure(list(fit = fit, G = G, N = N, Y = Y,
                       iter = stanargs$iter, chains = stanargs$chains,
                       thin = stanargs$thin, model_name = model_name), 
                  class = "bnlfa_fit")
  return(bm)
}

#' Extract the MAP pseudotime values from a \code{bnlfa_fit}
#' 
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' 
#' @export
#' 
#' @return MAP pseudotime vector of length N
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
  stopifnot(is(bm, "bnlfa_fit"))
  posterior.mode(mcmc(extract(bm$fit, "t")$t))
}

#' Reconstructed pseudotimes
#' @export
rexprs <- function(bm) UseMethod("rexprs")

#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
rexprs.bnlfa_fit <- function(bm) {
  stopifnot(is(bm, "bnlfa_fit"))
  Z <- apply(extract(bm$fit, "mu")$mu, 3, function(x) posterior.mode(mcmc(x)))
  Z <- t(Z)
  colnames(Z) <- colnames(bm$Y)
  return(Z)
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
#' Plot a \code{bnlfa_fit} object. Returns either a trace fit, MAP fit or MCMC diagnostic fit.
#' See the individual function calls (described below) for more details.
#' 
#' @param bm An object of class \code{bnlfa_fit}
#' @param what One of
#' \itemize{
#' \item \code{trace} This produces a heatmap of gene expression as a function of pseudotime
#' across different pseudotime samples. Underlying call is to \code{\link{plot_bnlfa_fit_trace}}.
#' \item \code{map} This plots gene expression as a function of the MAP pseudotime with a red
#' line denoting a LOESS fit (showing the overall trend). Underlying call is to
#' \code{\link{plot_bnlfa_fit_map}}
#' \item \code{diagnostic} This returns trace and autocorrelation plots of the log-posterior
#' probability. Underlying call is to \code{\link{plot_bnlfa_fit_diagnostics}}
#' }
#' @param ... Additional arguments passed to the corresponding functions
#' 
#' @return A \code{ggplot2} plot.
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
#' Plot MCMC diagnostics (traceplot and autocorrelation) of the log-posterior probability
#' for a \code{bnlfa_fit} object.
#' 
#' Further assessment of convergence can be done using \code{rstan} functions on the
#' underlying STAN object (accessed through \code{bm$fit}).
#' 
#' @param bm A \code{bnlfa_fit} object
#' @param nrow Number of rows. If 1, plots are side-by-side; if 2, plots are vertically aligned.
#' @export
#' @importFrom cowplot plot_grid
#' 
#' @return A \code{ggplot2} object
#' 
plot_bnlfa_fit_diagnostics <- function(bm, arrange = c("vertical", "horizontal")) {
  stopifnot(is(bm, "bnlfa_fit"))
  arrange <- match.arg(arrange)
  nrow <- switch(arrange,
                 vertical = 2,
                 horizontal = 1)
  plt <- cowplot::plot_grid(stan_trace(bm$fit, "lp__"), stan_ac(bm$fit, "lp__"), nrow = nrow)
  return(plt)
}

#' Plot heatmaps of gene expression changes
#' 
#' Produces a heatmap of gene expression as a function of pseudotime
#' across different pseudotime samples.
#' 
#' @param bm An object of class \code{bnlfa_fit}
#' @param samples Number of posterior pseudotime samples to use (number of rows of heatmap)
#' @param genes A vector that subsets the gene expression matrix. Defaults to the first \code{g}
#' genes, where \code{g} is either 4 or the number of genes in the model if less than 4.
#' @param output If \code{grid} then \code{cowplot::plot_grid} is called and a grid plot
#' of all genes is returned. If \code{plotlist} then a list of \code{ggplot2} objects is returned 
#' for the user to customise further
#' @param show_legend Logical. If \code{TRUE} then the legend (ie gene expression magnitude) is
#' displayed for each heatmap.
#' @param ... Additional arguments passed to \code{cowplot::plot_grid}
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
#' @return A \code{ggplot2} object.
plot_bnlfa_fit_trace <- function(bm, samples = 50, genes = seq_len(min(bm$G, 6)),
                                 output = c("grid", "plotlist"), 
                                 show_legend = FALSE, ...) {
  stopifnot(is(bm, "bnlfa_fit"))
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
      theme_bw() + 
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

#' Generic function to return sigmoid whenever needed
tsigmoid <- function(mu0, k, t0, t) {
  return( 2 * mu0 / (1 + exp(-k*(t - t0))))
}

#' Plot gene expression as a function of MAP pseudotime
#' 
#' Plot gene expression as a function of the MAP pseudotime with a red
#' line denoting a LOESS fit (showing the overall trend). Genes are plotted with
#' one per grid square (using a call to \code{facet_wrap(~ gene)}).
#' 
#' @param bm An object of class \code{bnlfa_fit}
#' @param genes A vector that subsets the gene expression matrix. Defaults to the first \code{g}
#' genes, where \code{g} is either 4 or the number of genes in the model if less than 4.
#' @param expression_units Units for expression to be displayed in along y-axis.
#' 
#' @importFrom reshape2 melt
#' @importFrom rstan extract
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda mcmc
#' @importFrom dplyr inner_join
#' @import ggplot2
#' 
#' @export
#' 
#' @return An object of class \code{ggplot2}
plot_bnlfa_fit_map <- function(bm, genes = seq_len(min(bm$G, 6)),
                               expression_units = "log2(TPM+1)") {
  stopifnot(is(bm, "bnlfa_fit"))
  tmap <- map_pseudotime(bm)
  Y <- bm$Y
  
  ## want to plot sigmoid function so need MAP estimates
  extr <- extract(bm$fit, pars = c("mu0", "k", "t0"))
  mu0_map <- posterior.mode(mcmc(extr$mu0))
  k_map <- posterior.mode(mcmc(extr$k))
  t0_map <- posterior.mode(mcmc(extr$t0))
  sig_map <- data.frame(mapply(tsigmoid, mu0_map, k_map, t0_map, MoreArgs = list(t = tmap)))
  names(sig_map) <- colnames(Y)
  
  ## Create data frame for gene expression values
  Y <- Y[,genes]
  dy <- data.frame(Y, pseudotime = tmap)
  dm <- melt(dy, id.vars = "pseudotime", 
             variable.name = "gene", 
             value.name = "expression")
  
  S <- sig_map[,genes]
  ds <- data.frame(S, pseudotime = tmap)
  dm2 <- melt(ds, id.vars = "pseudotime",
              variable.name = "gene",
              value.name = "predicted_expression")
  
  dm_joined <- inner_join(dm, dm2, by = c("pseudotime", "gene"))
  
  plt <- ggplot(dm_joined, aes(x = pseudotime, y = expression, colour = "Measured")) + 
    geom_point() +
    facet_wrap(~ gene, scales = "free_y") +
    xlab("MAP pseudotime") + ylab(paste("Expression", expression_units)) + theme_bw()
  plt <- plt + 
    geom_line(aes(x = pseudotime, y = predicted_expression, color = 'Predicted'), 
              size = 2, alpha = 0.7) +
    scale_colour_manual(values = c("Predicted" = "red", "Measured" = "black"), name = element_blank()) +
    theme(legend.position = "bottom")
  return( plt )
}

#' Plot heatmaps showing comparisons of measured data and imputed
#' 
#' @param return_plotlist If TRUE then the list of \code{ggplot}s is returned
#' instead of being plotted with \code{cowplot::plot_grid}
#' 
#' @export
#' @import ggplot2
#' 
#' @return Either a list of plots of class \code{ggplot} or a single 
#' \code{ggplot} showing them
plot_bnlfa_fit_comparison <- function(bm, return_plotlist = FALSE) {
  stopifnot(is(bm, "bnlfa_fit"))
  X <- bm$Y
  Z <- rexprs(bm)
  tmap <- map_pseudotime(bm)
  
  make_tile_plot <- function(X, tmap) {
    Xp <- apply(X, 2, function(x) (x - min(x)) / (max(x) - min(x)))
    
    dfx <- data.frame(Xp, pseudotime = rank(tmap)) %>%
      melt(id.vars = "pseudotime", variable.name = "gene", value.name = "expression")
    ggplot(dfx, aes(x = pseudotime, y = gene, fill = expression)) + 
      geom_tile() + scale_fill_viridis(name = expression("Relative\nexpression")) + 
      theme_bw() + 
      theme(axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title.y = element_blank()) +
      xlab("Pseudotime order")
  }
  plts <- lapply(list(X, Z), make_tile_plot, tmap)
  if(return_plotlist) return(plts)
  cowplot::plot_grid(plotlist = plts, nrow = 2, 
                     labels = c("Measured", "Reconstructed"))
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