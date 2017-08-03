# Ouija

Ouija is a statistical framework for learning interpretable pseudotimes from single-cell RNA-seq data using only small panels of marker genes. 

<img src="inst/www/fig_main.png" width="600"/>

Ouija uses nonlinear factor analysis (**A**) where the pseudotimes are represented by the latent variables. It models gene expression behaviour as either sigmoidal or transient (**B**) which in turn provides interpretable parameter estimates such as switch and peak times of genes. This allows us to formulate a Bayesian hypothesis test to work out whether a gene is regulated before another (**C**). We can also use the probablistic formulation to identify metastable states over pseudotime (**D**).

## Getting started

### Installation

```r
# install.packages("devtools")
devtools::install_github("kieranrcampbell/ouija")
```

To build the Ouija vignette install using

```r
devtools::install_github("kieranrcampbell/ouija", local = FALSE, 
                          args = "--preclean", build_vignettes = TRUE)
```

### Model fitting

Input is a cell-by-gene expression matrices that is non-negative and represents logged gene expression values. We recommend using `log2(TPM + 1)`. This can either take the form of a matrix or an `ExpressionSet` such as from the [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) package:

```r
library(ouija)
data(synth_gex) # synthetic gene expression data bundled
oui <- ouija(synth_gex)
pseudotimes <- map_pseudotime(oui)
```

For further usage options see the vignette.


## Authors

Kieran Campbell & Christopher Yau  
Wellcome Trust Centre for Human Genetics, University of Oxford

## Artwork

<img src="inst/www/chris_ouija.jpg" width="500"/>



