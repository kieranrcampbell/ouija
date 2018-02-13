
#' Fit descriptive single-cell pseudotime using Ouija
#' 
#' Given single-cell expression data fit a Bayesian non-linear factor analysis model using Ouija.
#' This infers pseudotimes along with interpretable gene behaviour parameters 
#' (corresponding to switch times and strengths, peak times and lengths). Each gene must be specified 
#' beforehand as switch-like or transient (default to all switch-like). Priors on all parameters
#' may also be specified.
#' 
#' \strong{Input format}
#' 
#' Ouija takes input in three formats:
#' 
#' \enumerate{
#' \item A cell-by-gene expression matrix of non-negative values. 
#' We recommend using log2(TPM + 1) or log2(RPKM + 1) as this 
#' is what the mean-variance relationship in the model is designed for.
#' \item A \code{SingleCellExperiment} (from the \pkg{SingleCellExperiment}) package (\strong{recommended})
#' \item An \code{ExpressionSet} (e.g. if you're using the \pkg{scater} package)
#' where \code{exprs(SCESet)} corresponds to the transpose of the matrix in (1) (kept for
#' backwards compatibility with older versions of \pkg{scater})
#' }
#' 
#' 
#' @param x Input expression. See details below.
#' @param switch_strengths Prior means of switch strengths
#' @param switch_times Prior means of switch times
#' @param switch_strength_sd Prior standard deviations of switch strengths
#' @param switch_time_sd Prior standard deviations of switch times
#' @param peak_times Prior means of peak times
#' @param peak_lengths Prior means of peak lengths
#' @param peak_time_sd Prior standard deviations of peak times
#' @param peak_length_sd Prior standard deviations of peak lengths
#' @param ... Additional arguments to \code{rstan::sampling} or \code{rstan::vb}
#' @param student_df Degrees of freedom for the student's t likelihood
#' @param inference_type The type of inference to be performed, either \code{hmc} for Hamiltonian
#' Monte Carlo or \code{vb} for ADVI (Variational Bayes). Note that HMC is typically more accurate
#' but VB will be orders of magnitude faster.
#' @param response_type A vector declaring whether each gene exhibits "switch" or "transient"
#' expression. Defaults to "switch" for all genes
#' @param single_cell_experiment_assay Character vector specifying the assay from 
#' \code{SingleCellExperiment} to use. 
#' Defaults to \code{"logcounts"}, so the input expression
#' matrix used is \code{assay(single_cell_experiment, "logcounts")}.
#' @param normalise_expression Logical, default TRUE. If TRUE the data is pre-normalised
#' so the average peak expression is approximately 1. This makes the strength parameters
#' approximately comparable between genes.
#' 
#' @import rstan
#' @importFrom Rcpp loadModule
#' @importFrom stats coef lm prcomp
#' 
#' @export
#' 
#' @return An object of type \code{ouija_fit} which contains both the \code{stan} fit
#' along with information about the setup. For more information see the vignette 
#' (online at \url{http://kieranrcampbell.github.io/ouija/}).
#' 
#' @examples 
#' \dontrun{
#' data(example_gex)
#' response_types <- c(rep("switch", 9), rep("transient", 2))
#' oui <- ouija(example_gex[1:40,], response_type = response_types, iter = 100)
#' }
ouija <- function(x,
                response_type = "switch",
                switch_strengths = NULL,
                switch_times = NULL,
                switch_strength_sd = NULL,
                switch_time_sd = NULL,
                peak_times = NULL,
                peak_lengths = NULL,
                peak_time_sd = NULL,
                peak_length_sd = NULL,
                student_df = 10,
                inference_type = c("hmc", "vb"),
                normalise_expression = TRUE,
                single_cell_experiment_assay = "logcounts",
                ...) {

  library(rstan)
  model_file <- "ouija.stan"

  inference_type <- match.arg(inference_type)

  Y <- NULL
  if (is(x, "ExpressionSet")) {
    ## convert to expression matrix Y
    Y <- t(Biobase::exprs(x))
  } else if (is(x, "SingleCellExperiment")) {
    Y <- t(x@assays[[ single_cell_experiment_assay ]])
  } else {
    Y <- x
  }
  if (!is(Y, "matrix")) {
    stop("x must either be an SCESet or matrix of gene expression values")
  }

  ## Now sanitize the input
  G <- ncol(Y) # number of genes
  N <- nrow(Y) # number of cells

  if (any(Y < 0)) {
    stop("Negative entries found in Y - Ouija supports non-negative expression")
  }

  if (student_df <= 0) {
    stop("Degrees of freedom of student distribution must be positive")
  }

  norm_factors <- rep(1, G)
  if (normalise_expression) {
    norm_factors <- apply(Y, 2, function(x) mean(x[x > 0]))
    Y <- t(t(Y) / norm_factors)
  }

  if (length(response_type) == 1) response_type <- rep(response_type, G)
  stopifnot(all(response_type %in% c("switch", "transient")))

  is_switch <- which(response_type == "switch")
  is_transient <- which(response_type == "transient")
  G_switch <- length(is_switch)
  G_transient <- length(is_transient)
  Y_switch <- Y[, is_switch]
  Y_transient <- Y[, is_transient]

  # we can fill in some values if they're null
  if (is.null(switch_strengths)) switch_strengths <- rep(0, G_switch)
  if (is.null(switch_strength_sd)) switch_strength_sd <- rep(5, G_switch)
  if (is.null(switch_times)) switch_times <- rep(0.5, G_switch)
  if (is.null(switch_time_sd)) switch_time_sd <- rep(1, G_switch)

  if (is.null(peak_times)) peak_times <- rep(0.5, G_transient)
  if (is.null(peak_time_sd)) peak_time_sd <- rep(0.1, G_transient)
  if (is.null(peak_lengths)) peak_lengths <- rep(50, G_transient)
  if (is.null(peak_length_sd)) peak_length_sd <- rep(10, G_transient)

  stopifnot(length(switch_strengths) == G_switch)
  stopifnot(length(switch_strength_sd) == G_switch)
  stopifnot(length(switch_times) == G_switch)
  stopifnot(length(switch_time_sd) == G_switch)

  stopifnot(length(peak_times) == G_transient)
  stopifnot(length(peak_time_sd) == G_transient)
  stopifnot(length(peak_lengths) == G_transient)
  stopifnot(length(peak_length_sd) == G_transient)

  ## stan setup
  data <- list(Y_switch = t(Y_switch), Y_transient = t(Y_transient),
               G_switch = G_switch, G_transient = G_transient, G = G, N = N,
               k_means = switch_strengths, k_sd = switch_strength_sd,
               t0_means = switch_times, t0_sd = switch_time_sd,
               p_means = peak_times, p_sd = peak_time_sd,
               b_means = peak_lengths, b_sd = peak_length_sd,
               student_df = student_df)


  stanfile <- system.file(model_file, package = "ouija")
  # model <- stan_model(stanfile, save_dso = TRUE, verbose = FALSE)

  ## manipulate stan defaults
  stanargs <- list(...)

  if (inference_type == "hmc") {
    if (!("iter" %in% names(stanargs))) stanargs$iter <- 1e4
    if (!("warmup" %in% names(stanargs))) stanargs$warmup <- stanargs$iter / 2
    if (!("chains" %in% names(stanargs))) stanargs$chains <- 1
    if (!("thin" %in% names(stanargs))) {
      # always specify thin so that approximately 1000 samples are returned
      stanargs$thin <- ceiling( ( stanargs$iter - stanargs$warmup ) / 1000)
    }
  }
  if (!("init" %in% names(stanargs))) {
    # We'll initialise to PC1 of switching genes if no init specified
    pc1 <- prcomp(Y_switch)$x[, 1]
    pc1_scaled <- (pc1 - min(pc1) + 1e-2) / (max(pc1) - min(pc1) + 2e-2)

    k_inits <- apply(Y_switch, 2, function(y) {
      coef(lm(y ~ pc1_scaled))[2]
    })

    if (inference_type == "hmc") {
      inits <- list()
      init_chain <- list(t = pc1_scaled, k = k_inits)
      for (i in seq_len(stanargs$chains)) inits[[i]] <- init_chain
      stanargs$init <- inits
    } else {
      stanargs$init <- list(t = pc1_scaled, k = k_inits)
    }
  }

  stanargs$data <- data


  ## call inference
  fit <- NULL
  if (inference_type == "hmc") {
    stanargs$file <- stanfile
    fit <- do.call(stan, stanargs)
  } else {
    model <- stan_model(stanfile, save_dso = TRUE, verbose = FALSE)
    stanargs$object <- model
    fit <- do.call(vb, stanargs)
  }

  oui <- structure(list(fit = fit,
                        G = G,
                        N = N,
                        Y = Y,
                        inference_type = inference_type,
                        iter = stanargs$iter,
                        chains = stanargs$chains,
                        thin = stanargs$thin,
                        strengths = switch_strengths,
                        switch_strength_sd = switch_strength_sd,
                        switch_times = switch_times,
                        switch_time_sd = switch_time_sd,
                        peak_times = peak_times,
                        peak_time_sd = peak_time_sd,
                        peak_lengths = peak_lengths,
                        peak_length_sd = peak_length_sd,
                        normalise_expression = normalise_expression,
                        norm_factors = norm_factors,
                        response_type = response_type),
                  class = "ouija_fit")
  return(oui)
}

#' Example gene expression matrix
#' 
#' A matrix containing some example gene expression data for 
#' 400 cells and 11 genes, the first 9 of which exhibit switch-like
#' expression and the 
#' 
#' @return A 400-by-11 example expression matrix
#' 
#' @examples 
#' data(example_gex)
"example_gex"


#' Precomputed Ouija fit
#' 
#' An example \code{ouija_fit} to avoid
#' large times computing while testing.
#' 
#' @seealso example_gex
#' @examples
#' data(oui)
"oui"
