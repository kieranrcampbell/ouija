# BNLFA

R package for Bayesian non-linear factor analysis for incorporating prior knowledge in single-cell trajectory learning.

## Warning

If there was a greek letter before alpha, it would describe this package. 

## Installation

```R
# install.packages("devtools")
devtools::install_github("kieranrcampbell/bnlfa")
```

## Usage

```R
## Using gene expression matrix Y
bm <- bnlfa(Y, prior_info, model_type, mcmc_pars)

## Using scater object
bm <- bnlfa(sce, prior_info, model_type, mcmc_pars)
```

## Authors

Kieran Campbell & Christopher Yau  
Wellcome Trust Centre for Human Genetics, University of Oxford

