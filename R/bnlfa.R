## Bayesian non-linear factor analysis for incorporating 
## prior knowledge into single-cell trajectory learning
## kieranc@well.ox.ac.uk

library(MCMCglmm)
library(coda)
library(rstan)



#' Fit Bayesian non-linear factor analysis
#' 
#' @param Y Cell-by-gene (N by G) logged expression matrix
#' @param k_means G mean activation strength parameters
#' @param t0_means G mean activation time parameters
#' @param k_sd Optional standard deviations for k parameters
#' @param t0_sd Optional standard deviations for t0 parameters
#' 
#' @returns ...
bnlfa <- function(Y, k_means, t0_means,
                  k_sd = rep(0.5, ncol(Y)), 
                  t0_sd = rep(0.5, ncol(Y)),
                  model_name = "bnlfa.stan",
                  ...) {
  G <- ncol(Y)
  N <- nrow(Y)
  
  ## sanity checking
  stopifnot(length(k_means) == length(t0_means))
  stopifnot(length(k_means) == G)
  
  ## stan setup
  data <- list(Y = t(Y), G = G, N = N,
               k_means = k_means, k_sd = k_sd,
               t0_means = t0_means, t0_sd = t0_sd)
  
  stanfile <- system.file(model_name, package = "bnlfa")
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
  
  ## get pseudotime map
  tmap <- posterior.mode(mcmc(extract(fit, "t")$t))
  return(list(fit = fit, tmap = tmap))
}

#' Print bnfla model
print_model <- function() {
  stanfile <- system.file("bnlfa.stan", package = "bnlfa")
  model <- stan_model(stanfile)
  print(model)
}