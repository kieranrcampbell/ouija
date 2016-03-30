
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user

}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0[G];
  real<lower = 0> tau[G]; // precision
  
  // parameters with user-defined priors
  real k[G];
  real t0[G];
  
  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell

  // hyper-parameters for variances
  real<lower = 0> alpha;
  real<lower = 0> beta;
  
  real nu; // hyper-parameter for means
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

  nu ~ exponential(1);
  mu0 ~ normal(nu, 1);
  alpha ~ gamma(2,1);
  beta ~ gamma(2,1);
  
  for(g in 1:G) tau[g] ~ gamma(alpha, beta);
  for(i in 1:N) t[i] ~ normal(0.5, 1);
  
  for(g in 1:G) {
    for(i in 1:N) {
      Y[g][i] ~ normal(mu[g][i], 1 / sqrt(tau[g]));
    }
  }
}