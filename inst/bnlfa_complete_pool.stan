
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user

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
  for(g in 1:G) {
    for(i in 1:N) mu[g][i] <- 2 * mu0[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
  }
}

model {
  // user defined priors
  for(g in 1:G) {
    k[g] ~ normal(k_means[g], k_sd[g]);
    t0[g] ~ normal(t0_means[g], t0_sd[g]);
  }
  
  // model priors

  mu0 ~ exponential(1);
  tau ~ gamma(2.0, 1.0);
  for(i in 1:N) t[i] ~ normal(0.5, 1);
  
  for(g in 1:G) {
    for(i in 1:N) {
      Y[g][i] ~ normal(mu[g][i], 1 / sqrt(tau));
    }
  }
}