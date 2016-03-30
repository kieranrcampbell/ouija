
/* Bayesian nonlinear factor analysis using gene directions only
kieran.campbell@sjc.ox.ac.uk */

data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  // direction parameters (ie k = -1 or k = 1) for all genes
  int<lower = -1, upper = 1> k_bit[G]; 
}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0[G];
  real<lower = 0> tau; // precision
  
  // parameters with user-defined priors
  real<lower = 0> k[G];
  // real<lower = 0> kdelt[G]; // ARD for each gene
  real t0[G];
  
  real<lower = 0> k_prior;
  real<lower = 0> mu0_prior;

  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell
}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  for(g in 1:G) {
    for(i in 1:N) mu[g][i] <- 2 * mu0[g] / (1 + exp(-k_bit[g] * k[g] * (t[i] - t0[g])));
  }
}

model {
  // user defined priors
  for(g in 1:G) {
    k[g] ~ exponential(k_prior);
    mu0[g] ~ exponential(mu0_prior);
    t0[g] ~ normal(0, 1);
  }
  
  // model priors
  k_prior ~ gamma(2,1);
  mu0_prior ~ gamma(2,1);

  // for(g in 1:G) tau[g] ~ gamma(2.0, 1.0);
  tau ~ gamma(2.0, 1.0);
  for(i in 1:N) t[i] ~ normal(0.5, 1);
  
  for(g in 1:G) {
    for(i in 1:N) {
      Y[g][i] ~ normal(mu[g][i], 1 / sqrt(tau));
    }
  }
}