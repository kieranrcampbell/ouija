
data {
  int<lower = 2> N; // number of cells
  int<lower = 2> G; // number of genes
  
  vector<lower = 0>[N] Y[G]; // matrix of gene expression values
  
  real k_means[G]; // mean parameters for k provided by user
  real k_sd[G]; // standard deviation parameters for k provided by user
  
  real t0_means[G]; // mean parameters for t0 provided by user
  real t0_sd[G]; // standard deviation parameters for t0 provided by user
  
  real student_df;
}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0[G];
  
  real<lower = 0> phi; // mean-variance "overdispersion" parameter 
  
  // parameters with user-defined priors
  real k[G];
  real t0[G];
  
  real<lower = 0> mu_hyper;
  
  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell

  real beta[2];
}

transformed parameters {
  vector[N] mu[G]; // mean for cell i gene g
  vector<lower = 0>[N] ysd[G];

  for(g in 1:G) {
    for(i in 1:N) {
      mu[g][i] = 2 * mu0[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
      ysd[g][i] = sqrt( (1 + phi) * mu[g][i] + 0.01);
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
  phi ~ gamma(12, 4); 
  
  t ~ normal(0.5, 1);

  beta ~ normal(0, 0.1);
  
  // Zero inflation per-mean
  for(g in 1:G) {
    for(i in 1:N) {
      if(Y[g][i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | beta[1] + beta[2] * mu[g][i]),
                              bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu[g][i]) + 
                              student_t_lpdf(Y[g][i] | student_df, mu[g][i], ysd[g][i]));
      } else {
        target += bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu[g][i]) + 
        student_t_lpdf(Y[g][i] | student_df, mu[g][i], ysd[g][i]);
      }
    }
  }
}

