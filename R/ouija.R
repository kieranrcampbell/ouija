## Bayesian non-linear factor analysis for incorporating 
## prior knowledge into single-cell trajectory learning
## kieranc@well.ox.ac.uk


#' Fit a Ouija object.
#' 
#' Fit a Bayesian non-linear factor analysis model given some single-cell
#' gene expression data.
#' 
#' This function takes either a \code{SCESet} object or an expression matrix
#' and returns a \code{ouija_fit} object including posterior traces for all
#' variables.
#' 
#' @param x Either an \code{SCESet} from \code{scater} or a
#' cell-by-gene (N by G) matrix of non-negative values representing gene expression.
#' log2(TPM + 1) is recommended.
#' @param strengths G mean activation strength parameters
#' @param times G mean activation time parameters
#' @param strength_sd Optional standard deviations for k parameters
#' @param time_sd Optional standard deviations for t0 parameters
#' @param warn_lp Ouija can perform a crude check of convergence in cases where may
#' models are being fit and manual inspection may be cumbersome. The log-likelihood after
#' the burn period is regressed off the iteration number, and if the gradient of the fit
#' falls above a threshold (set by \code{lp_gradient_threshold}) then the user is warned.
#' @param lp_gradient_threshold The threshold for convergence warning. If the slope of regressing
#' the log-probability of the model against the iteration number falls above this value then
#' the user is warned.
#' @param ... Additional arguments to \code{rstan::sampling}
#' @param student_df Degrees of freedom for the student's t likelihood
#' @param inference_type The type of inference to be performed, either \code{hmc} for Hamiltonian
#' Monte Carlo or \code{vb} for ADVI (Variational Bayes). Note that HMC is typically more accurate
#' but VB will be orders of magnitude faster.
#' 
#' @param normalise_expression Logical, default TRUE. If TRUE the data is pre-normalised
#' so the average peak expression is approximately 1. This makes the strength parameters
#' approximately comparable between genes.
#' 
#' @import rstan
#' @import stats
#' @importFrom Rcpp loadModule
#' 
#' @export
#' 
#' @return An object of type \code{ouija_fit}
#' 
#' @examples 
#' data(synth_gex)
#' oui <- ouija(synth_gex, strengths = 5 * c(1, -1, 1, -1, -1, -1), iter = 100)
ouija <- function(x, 
                  response_type = "switch",
                  strengths = NULL, times = NULL,
                  strength_sd = NULL, time_sd = NULL,
                  peak_times = NULL, bandwidths = NULL,
                  peak_sd = NULL, bandwidth_sd = NULL,
                  warn_lp = TRUE,
                  lp_gradient_threshold = 1e-2,
                  student_df = 10,
                  inference_type = c("hmc", "vb"),
                  normalise_expression = TRUE,
                  ...) {
  
  # requireNamespace('rstan')
  model_file <- "ouija.stan"
  
  inference_type <- match.arg(inference_type)

  Y <- NULL
  if(is(x, "SCESet")) {
    ## convert to expression matrix Y  
    Y <- t(Biobase::exprs(x))
  } else {
    Y <- x
  }
  if(!is(Y, "matrix")) {
    stop("x must either be an SCESet or matrix of gene expression values")
  }
  
  ## Now sanitize the input
  G <- ncol(Y) # number of genes
  N <- nrow(Y) # number of cells
  
  if(any(Y < 0)) {
    stop("Negative entries found in Y - Ouija supports non-negative expression")
  }
  
  if(student_df <= 0) {
    stop("Degrees of freedom of student distribution must be positive")
  }
  
  ## -- Normalise the dataset -- ##
  if(normalise_expression) {
    norm_factors <- apply(Y, 2, function(x) mean(x[x > 0]))
    Y <- t(t(Y) / norm_factors)
  }
  
  ## Sort out switch-vs-transient stuff
  if(length(response_type) == 1) response_type <- rep(response_type, G)
  stopifnot(all(response_type %in% c("switch", "transient")))
  
  is_switch <- which(response_type == "switch")
  is_transient <- which(response_type == "transient")
  G_switch <- length(is_switch)
  G_transient <- length(is_transient)
  Y_switch <- Y[,is_switch]
  Y_transient <- Y[,is_transient]
  
  # we can fill in some values if they're null
  if(is.null(strengths)) strengths = rep(0, G_switch)
  if(is.null(strength_sd)) strength_sd <- rep(1, G_switch)
  if(is.null(times)) times <- rep(0.5, G_switch) ## change if constrained
  if(is.null(time_sd)) time_sd <- rep(1, G_switch)
  
  if(is.null(peak_times)) peak_times <- rep(0.5, G_transient)
  if(is.null(peak_sd)) peak_sd <- rep(0.1, G_transient)
  if(is.null(bandwidths)) bandwidths <- rep(50, G_transient)
  if(is.null(bandwidth_sd)) bandwidth_sd <- rep(10, G_transient)

  stopifnot(length(strengths) == G_switch)
  stopifnot(length(strength_sd) == G_switch)
  stopifnot(length(times) == G_switch)
  stopifnot(length(time_sd) == G_switch)
  
  stopifnot(length(peak_times) == G_transient)
  stopifnot(length(peak_sd) == G_transient)
  stopifnot(length(bandwidths) == G_transient)
  stopifnot(length(bandwidth_sd) == G_transient)

  ## stan setup
  data <- list(Y_switch = t(Y_switch), Y_transient = t(Y_transient),
               G_switch = G_switch, G_transient = G_transient, G = G, N = N,
               k_means = strengths, k_sd = strength_sd,
               t0_means = times, t0_sd = time_sd,
               p_means = peak_times, p_sd = peak_sd,
               b_means = bandwidths, b_sd = bandwidth_sd,
               student_df = student_df)
  
  
  stanfile <- system.file(model_file, package = "ouija")
  model <- stan_model(stanfile, save_dso = FALSE, verbose = TRUE)
  
  ## manipulate stan defaults
  stanargs <- list(...)
  if(inference_type == "hmc") { # These are all MCMC parameters
    if(!('iter' %in% names(stanargs))) stanargs$iter <- 1e4
    if(!('warmup' %in% names(stanargs))) stanargs$warmup <- stanargs$iter / 2
    if(!('chains' %in% names(stanargs))) stanargs$chains <- 1
    if(!('thin' %in% names(stanargs))) {
      # always specify thin so that approximately 1000 samples are returned
      stanargs$thin <- ceiling((stanargs$iter - stanargs$warmup) / 1000)
    }
    if(!('init' %in% names(stanargs))) {
      ## We'll initialise to PC1 of switching genes if no init specified
      pc1 <- prcomp(Y_switch)$x[,1]
      pc1_scaled <- (pc1 - min(pc1) + 1e-2) / (max(pc1) - min(pc1) + 2e-2)
      init_chain <- list(t = pc1_scaled)
      inits <- list()
      for(i in seq_len(stanargs$chains)) inits[[i]] <- init_chain
      stanargs$init <- inits
    }
  }
  stanargs$object <- model
  stanargs$data <- data
  
  ## call inference
  fit <- NULL
  if(inference_type == "hmc") {
    fit <- do.call(sampling, stanargs)
  } else {
    fit <- do.call(vb, stanargs)
  }
  
  
  ## Do a really dumb automated check of convergence:
  # if(warn_lp && inference_type == "hmc") {
  #   lp <- extract(fit, pars = "lp__")$lp__
  #   siter <- seq_along(lp)
  #   lplm <- lm(lp ~ siter)
  #   if(coef(lplm)[2] > lp_gradient_threshold) {
  #     warning(paste("Gradient of log-probability against iteration greater than threshold: "), coef(lplm)[2])
  #     warning("Model may not be converged")
  #   }
  # }
    
  oui <- structure(list(fit = fit, G = G, N = N, Y = Y,
                        inference_type = inference_type,
                        iter = stanargs$iter, chains = stanargs$chains,
                        thin = stanargs$thin,
                        strengths = strengths, strength_sd = strength_sd,
                        times = times, time_sd = time_sd,
                        normalise_expression = normalise_expression,
                        norm_factors = norm_factors,
                        response_type = response_type), 
                  class = "ouija_fit")
  return(oui)
}

#' Synthetic gene expression matrix
#' 
#' A matrix containing some synthetic gene expression data for 
#' 100 cells and 6 genes
#' 
#' @return A matrix containing some synthetic gene expression data for 
#' 100 cells and 6 genes
#' 
#' @examples 
#' data(synth_gex)
"synth_gex"

#' Synthetic gene pseudotime vector
#' 
#' A vector with the 'true' pseudotimes for the synthetic 
#' gene expression data in \code{synth_gex}
#' 
#' @return A vector with the 'true' pseudotimes for the synthetic 
#' gene expression data in \code{synth_gex}
#' 
#' @examples
#' data(true_pst)
"true_pst"

#' Precomputed Ouija fit
#' 
#' The result of calling \code{oui <- ouija(synth_gex)} to avoid
#' large times computing while testing.
#' 
#' @seealso synth_gex
#' @examples
#' data(oui)
"oui"
