
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user

  real<lower = 0> lambda; // hyper-hyper parameter for nu
}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0[G];
  real<lower = 0> tau[G]; // precision
  
  
  // parameters with user-defined priors
  real k[G];
  real t0[G];
  
  real<lower = 0> nu;
  
  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell


}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  real<lower = 0> one_over_sqrt_tau[G];
  vector<lower = 0>[N] ysd[G];

  for(g in 1:G) {
    one_over_sqrt_tau[g] <- 1 / sqrt(tau[g]);
    for(i in 1:N) {
      mu[g][i] <- 2 * mu0[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
      ysd[g][i] <- one_over_sqrt_tau[g] * (0.01 + mu[g][i]);
    }
    
  }
  
}

model {
  // user defined priors
  k ~ normal(k_means, k_sd);
  t0 ~ normal(t0_means, t0_sd);
  
  // model priors
  mu0 ~ exponential(1);
  tau ~ gamma(nu / 2, 2);
  nu ~ exponential(lambda);
  t ~ normal(0.5, 1);
  
  for(g in 1:G) {
    Y[g] ~ normal(mu[g], ysd[g]);
  }
}