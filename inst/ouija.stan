
data {
  int<lower = 2> N; // number of cells
  int<lower = 0> G_switch; // number of switch-like genes
  int<lower = 0> G_transient; // number of transient genes
  int<lower = 1> G;
  
  vector<lower = 0>[N] Y_switch[G_switch]; // matrix of gene expression values
  vector<lower = 0>[N] Y_transient[G_transient]; // matrix of gene expression values
  
  real k_means[G_switch]; // mean parameters for k provided by user
  real<lower = 0> k_sd[G_switch]; // standard deviation parameters for k provided by user
  
  real t0_means[G_switch]; // mean parameters for t0 provided by user
  real<lower = 0> t0_sd[G_switch]; // standard deviation parameters for t0 provided by user
  
  real p_means[G_transient]; // mean peak time
  real<lower = 0> p_sd[G_transient]; // sd peak time
  
  real b_means[G_transient]; // mean bandwidth
  real<lower = 0> b_sd[G_transient]; // sd bandwidth
  
  real student_df; // student d-o-f
}

parameters {
  // parameters we'll let stan infer
  real<lower = 0> mu0_switch[G_switch];
  real<lower = 0> mu0_transient[G_transient];
  
  real<lower = 0> phi; // mean-variance "overdispersion" parameter 
  
  // switch parameters
  real k[G_switch];
  real<lower = 0, upper = 1> t0[G_switch];
  
  // transient parameters
  real<lower = 0, upper = 1> p[G_transient];
  real<lower = 0> b[G_transient];
  
  real<lower = 0> mu_hyper;
  
  real<lower = 0, upper = 1> t[N]; // pseudotime of each cell

  real beta[2];
}

transformed parameters {
  vector[N] mu_switch[G_switch]; 
  vector<lower = 0>[N] sd_switch[G_switch];
  vector[N] mu_transient[G_transient]; 
  vector<lower = 0>[N] sd_transient[G_transient];


  for(g in 1:G_switch) {
    for(i in 1:N) {
      mu_switch[g][i] = 2 * mu0_switch[g] / (1 + exp(-k[g] * (t[i] - t0[g])));
      sd_switch[g][i] = sqrt( (1 + phi) * mu_switch[g][i] + 0.01);
    }
  }
  
  for(g in 1:G_transient) {
    for(i in 1:N) {
      mu_transient[g][i] = 2 * mu0_transient[g] * exp(- 10 * b[g] * (t[i] - p[g])^2);
      sd_transient[g][i] = sqrt( (1 + phi) * mu_transient[g][i] + 0.01);
    }
  }
}

model {
  // user defined priors
  k ~ normal(k_means, k_sd);
  t0 ~ normal(t0_means, t0_sd);
  
  p ~ normal(p_means, p_sd);
  b ~ normal(b_means, b_sd);
  
  // model priors
  mu0_switch ~ gamma(mu_hyper / 2, 0.5);
  mu0_transient ~ gamma(mu_hyper / 2, 0.5);
  
  /* Try tau ~ chi-sq(nu) using gamma approximation. See
  http://www.stat.umn.edu/geyer/old03/5102/notes/brand.pdf */
  phi ~ gamma(12, 4); 
  
  t ~ normal(0.5, 1);

  beta ~ normal(0, 0.1);
  
  // Switch likelihood
  for(g in 1:G_switch) {
    for(i in 1:N) {
      if(Y_switch[g][i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | beta[1] + beta[2] * mu_switch[g][i]),
                              bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu_switch[g][i]) + 
                              student_t_lpdf(Y_switch[g][i] | student_df, mu_switch[g][i], sd_switch[g][i]));
      } else {
        target += bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu_switch[g][i]) + 
        student_t_lpdf(Y_switch[g][i] | student_df, mu_switch[g][i], sd_switch[g][i]);
      }
    }
  }

  // Transient likelihood
  for(g in 1:G_transient) {
    for(i in 1:N) {
      if(Y_transient[g][i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | beta[1] + beta[2] * mu_transient[g][i]),
                              bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu_transient[g][i]) + 
                              student_t_lpdf(Y_transient[g][i] | student_df, mu_transient[g][i], sd_transient[g][i]));
      } else {
        target += bernoulli_logit_lpmf(0 | beta[1] + beta[2] * mu_transient[g][i]) + 
        student_t_lpdf(Y_transient[g][i] | student_df, mu_transient[g][i], sd_transient[g][i]);
      }
    }
  }
}

