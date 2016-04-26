
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user

  real<lower = 0> lambda; // hyper-hyper parameter for nu
  
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
  
  
  // parameters with user-defined priors
  real k[G];
  real t0[G];
  
  real<lower = 0> nu;
  real<lower = 0> mu_hyper;
  
  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell
  //real t[N];
  
  // real<lower = 0, upper = 1> theta[G];
  //real<lower = 0, upper = 1> theta[N];
  real beta[2];
}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  real<lower = 0> one_over_sqrt_tau[G];
  vector<lower = 0>[N] ysd[G];

  for(g in 1:G) {
    one_over_sqrt_tau[g] <- 1 / sqrt(tau[g]);
    for(i in 1:N) {
      mu[g][i] <- 2 * mu0[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
      if(mean_variance == 1) {
       # ysd[g][i] <- sqrt(mu[g][i] / tau[g] + 1e-4);
        ysd[g][i] <- sqrt(0.01 + mu[g][i] * (1 + mu[g][i] / tau[g]));
      } else {
        ysd[g][i] <- one_over_sqrt_tau[g];
      }
    }
  }
}

model {
  // user defined priors
  k ~ normal(k_means, k_sd);
  t0 ~ normal(t0_means, t0_sd);
  
  // model priors
  mu0 ~ gamma(mu_hyper / 2, 0.5);
  
  /* Try tau ~ chi-sq(nu) using gamma approximation. See
  http://www.stat.umn.edu/geyer/old03/5102/notes/brand.pdf */
  tau ~ gamma(nu / 2, 0.5);
  nu ~ exponential(lambda);
  //t ~ normal(0.5, 1);
  t ~ normal(0.5, 1);
  
  // for(i in 1:N) theta[i] ~ beta(2, 2);
  
  // // Zero inflation per-gene
  // for(g in 1:G) {
  //   for(i in 1:N) {
  //     if(Y[g][i] == 0) {
  //       increment_log_prob(log_sum_exp(bernoulli_log(1, theta[i]),
  //       bernoulli_log(0, theta[i]) + student_t_log(Y[g][i], 2, mu[g][i], ysd[g][i])));
  //     } else {
  //       increment_log_prob(bernoulli_log(0, theta[i]) +
  //       student_t_log(Y[g][i], 2, mu[g][i], ysd[g][i]));
  //     }
  //   }
  // }
  
  beta ~ normal(0, 0.1);
  
  // Zero inflation per-mean
  for(g in 1:G) {
    for(i in 1:N) {
      if(Y[g][i] == 0) {
        increment_log_prob(log_sum_exp(bernoulli_logit_log(1, beta[1] + beta[2] * mu[g][i]),
        bernoulli_logit_log(0, beta[1] + beta[2] * mu[g][i]) + 
        student_t_log(Y[g][i], 2, mu[g][i], ysd[g][i])));
      } else {
        increment_log_prob(bernoulli_logit_log(0, beta[1] + beta[2] * mu[g][i]) + 
        student_t_log(Y[g][i], 2, mu[g][i], ysd[g][i]));
      }
    }
  }
}