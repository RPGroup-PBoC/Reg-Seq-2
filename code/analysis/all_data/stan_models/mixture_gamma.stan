data {
  int<lower=0> N;
  array[N] real<lower=0> n;
  real alpha_mu; 
  real alpha_sigma;
  real beta_mu; 
  real beta_sigma;
}


parameters {
  vector[2] log10_alpha;
  vector[2] log10_beta;
  real<lower=0, upper=1> w;
}


transformed parameters {
  vector[2] alpha = 10 .^ log10_alpha;
  vector[2] beta_ = 10 .^ log10_beta;
}


model {
  // Priors
  log10_alpha ~ normal(alpha_mu, alpha_sigma);
  log10_beta ~ normal(beta_mu, beta_sigma);
  w ~ beta(1.0, 1.0);

  // Likelihood
  for (n_val in n) {
    target += log_mix(
      w,
      gamma_lupdf(n_val | alpha[1], beta_[1]),
      gamma_lupdf(n_val | alpha[2], beta_[2])
    );
  }
}


generated quantities {
  array[N] real n_ppc;
  array[N] real log_lik;

  // Posterior predictive checks
  for (i in 1:N) {
    if (uniform_rng(0.0, 1.0) < w) {
      n_ppc[i] = gamma_rng(alpha[1], beta_[1]);
    }
    else {
      n_ppc[i] = gamma_rng(alpha[2], beta_[2]);
    }
  }

  // Compute pointwise log likelihood
  for (i in 1:N) {
    log_lik[i] = log_mix(
      w,
      gamma_lpdf(n[i] | alpha[1], beta_[1]),
      gamma_lpdf(n[i] | alpha[2], beta_[2]));
    
  }
}