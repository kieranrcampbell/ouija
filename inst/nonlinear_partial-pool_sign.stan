
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values


  real<lower = 0> lambda; // hyper-hyper parameter for nu
  
  // direction parameters (ie k = -1 or k = 1) for all genes
  int<lower = -1, upper = 1> sign_bits[G];
  
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
  real<lower = 0> tau[G]; // precision
  
  real<lower = 0> k[G];
  real<lower = 0> k_prior;
  real t0[G];
  
  real<lower = 0> nu;
  
  # real<lower = 0, upper = 1> t[N]; // pseudotime of each cell
  real t[N];
}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  real<lower = 0> one_over_sqrt_tau[G];
  vector<lower = 0>[N] ysd[G];

  for(g in 1:G) {
    one_over_sqrt_tau[g] <- 1 / sqrt(tau[g]);
    for(i in 1:N) {
      mu[g][i] <- 2 * mu0[g] / (1 + exp(-sign_bits[g] * k[g] * (t[i] - t0[g])));
      if(mean_variance == 1) {
        ysd[g][i] <- sqrt(0.1 + mu[g][i] * (1 + mu[g][i] / tau[g]));
      } else {
        ysd[g][i] <- one_over_sqrt_tau[g];
      }
    }
    
  }
  
}

model {
  // user defined priors
  k ~ exponential(k_prior);
  t0 ~ normal(0.5, 1);
  k_prior ~ gamma(2,1);
  
  
  // model priors
  mu0 ~ exponential(1);
  tau ~ gamma(nu / 2, 2);
  nu ~ exponential(lambda);
  t ~ normal(0, 1);
  
  for(g in 1:G) {
    // Y[g] ~ student_t(1, mu[g], ysd[g]);
    Y[g] ~ normal(mu[g], ysd[g]);
  }
}