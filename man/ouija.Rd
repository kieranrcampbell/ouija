% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ouija.R
\name{ouija}
\alias{ouija}
\title{Fit descriptive single-cell pseudotime using Ouija}
\usage{
ouija(x, response_type = "switch", switch_strengths = NULL,
  switch_times = NULL, switch_strength_sd = NULL,
  switch_time_sd = NULL, peak_times = NULL, peak_lengths = NULL,
  peak_time_sd = NULL, peak_length_sd = NULL, student_df = 10,
  inference_type = c("hmc", "vb"), normalise_expression = TRUE,
  single_cell_experiment_assay = "logcounts", ...)
}
\arguments{
\item{x}{Input expression. See details below.}

\item{response_type}{A vector declaring whether each gene exhibits "switch" or "transient"
expression. Defaults to "switch" for all genes}

\item{switch_strengths}{Prior means of switch strengths}

\item{switch_times}{Prior means of switch times}

\item{switch_strength_sd}{Prior standard deviations of switch strengths}

\item{switch_time_sd}{Prior standard deviations of switch times}

\item{peak_times}{Prior means of peak times}

\item{peak_lengths}{Prior means of peak lengths}

\item{peak_time_sd}{Prior standard deviations of peak times}

\item{peak_length_sd}{Prior standard deviations of peak lengths}

\item{student_df}{Degrees of freedom for the student's t likelihood}

\item{inference_type}{The type of inference to be performed, either \code{hmc} for Hamiltonian
Monte Carlo or \code{vb} for ADVI (Variational Bayes). Note that HMC is typically more accurate
but VB will be orders of magnitude faster.}

\item{normalise_expression}{Logical, default TRUE. If TRUE the data is pre-normalised
so the average peak expression is approximately 1. This makes the strength parameters
approximately comparable between genes.}

\item{single_cell_experiment_assay}{Character vector specifying the assay from 
\code{SingleCellExperiment} to use. 
Defaults to \code{"logcounts"}, so the input expression
matrix used is \code{assay(single_cell_experiment, "logcounts")}.}

\item{...}{Additional arguments to \code{rstan::sampling} or \code{rstan::vb}}
}
\value{
An object of type \code{ouija_fit} which contains both the \code{stan} fit
along with information about the setup. For more information see the vignette 
(online at \url{http://kieranrcampbell.github.io/ouija/}).
}
\description{
Given single-cell expression data fit a Bayesian non-linear factor analysis model using Ouija.
This infers pseudotimes along with interpretable gene behaviour parameters 
(corresponding to switch times and strengths, peak times and lengths). Each gene must be specified 
beforehand as switch-like or transient (default to all switch-like). Priors on all parameters
may also be specified.
}
\details{
\strong{Input format}

Ouija takes input in three formats:

\enumerate{
\item A cell-by-gene expression matrix of non-negative values. 
We recommend using log2(TPM + 1) or log2(RPKM + 1) as this 
is what the mean-variance relationship in the model is designed for.
\item [\strong{Recommended}] A \code{SingleCellExperiment} (from the \pkg{SingleCellExperiment}) package 
}

\strong{Controlling inference with stan}

The \code{...} argument accepts various parameters than can be passed to \pkg{rstan} 
to control inference using HMC sampling. Specifically these are
\enumerate{
 \item \code{iter} Number of HMC iterations (default \code{1e4})
 \item \code{warmup} Number of iterations to use as warmup (defaults to \code{iter/2})
 \item \code{chains} Number of HMC chains to run in parallel (defaults to 1, but more recommended)
 \item \code{thin} The thinning of the markov chain (ie every \code{thin} iteration is kept). Defaults 
 so that 1000 samples are returned.
}
}
\examples{
\dontrun{
data(example_gex)
response_types <- c(rep("switch", 9), rep("transient", 2))
oui <- ouija(example_gex[1:40,], response_type = response_types, iter = 100)
}
}
