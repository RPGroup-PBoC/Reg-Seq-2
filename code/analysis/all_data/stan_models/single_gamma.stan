data {
  int<lower=0> N;
  array[N] real<lower=0> n;
  real alpha_mu; 
  real alpha_sigma;
  real beta_mu; 
  real beta_sigma;
}


parameters {
  real log10_alpha;
  real log10_beta;
}


transformed parameters {
  real alpha = 10^log10_alpha;
  real beta_ = 10^log10_beta;
}


model {
  // Priors
  log10_alpha ~ normal(alpha_mu, alpha_sigma);
  log10_beta ~ normal(beta_mu, beta_sigma);

  // Likelihood
  n ~ gamma(alpha, beta_);
}

generated quantities {
  array[N] real n_ppc;
  array[N] real log_lik;

  // Posterior predictive checks
  for (i in 1:N) {
    n_ppc[i] = gamma_rng(alpha, beta_);
  }

  // Compute pointwise log likelihood
  for (i in 1:N) {
    log_lik[i] = gamma_lpdf(n[i] | alpha, beta_);
  }
}