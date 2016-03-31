
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user

  /*
  Binary indicator variable to model the mean variance relationship:
  if 1 then variance = mean / tau
  if 0 then variance = 1 / tau
  */
  int<lower = 0, upper = 1> mean_variance;
}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0[G];
  real<lower = 0> tau; // precision
  
  // parameters with user-defined priors
  real k[G];
  real t0[G];

  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell
}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  real<lower = 0> one_over_sqrt_tau;
  vector<lower = 0>[N] ysd[G];

  for(g in 1:G) {
    for(i in 1:N) mu[g][i] <- 2 * mu0[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
  }
  
  one_over_sqrt_tau <- 1 / sqrt(tau);
  if(mean_variance == 1) {
    for(g in 1:G) ysd[g] <- one_over_sqrt_tau * (0.01 + mu[g]);
  } else {
    for(g in 1:G) ysd[g] <- one_over_sqrt_tau * rep_vector(1.0, N);
  }
}

model {
  // user defined priors
  k ~ normal(k_means, k_sd);
  t0 ~ normal(t0_means, t0_sd);
  
  // model priors
  mu0 ~ exponential(1);
  tau ~ gamma(2, 1);
  t ~ normal(0.5, 1);
  
  for(g in 1:G) {
    Y[g] ~ normal(mu[g], ysd[g]);
  }
}