data {
  int<lower=1> Ntotal;             // Total number of patients (rows in design matrix)
  int<lower=1> Ncol;               // Total number of columns in the design matrix
  matrix[Ntotal, Ncol] X;          // Design matrix (Doses 1-5 one-hot encoded + Covariates)
  array[Ntotal] int<lower=0, upper=1> y;  // DLT indicator (0 or 1)
  vector[Ntotal] w;                // TiTE observational weights (t_i / T_max)
  vector[5] prior_means;           // Logit-skeleton anchors for the 5 doses
}

parameters {
  vector[Ncol] betas;              // Combined coefficient vector: [alpha1..5, beta_pk, beta_age...]
}

model {
  // 1. Doses remain locked to your clinical BWS skeleton
  for (d in 1:5) {
    betas[d] ~ normal(prior_means[d], 0.2);
  }
  
  // 2. ALL Covariates (Categorical or Continuous) get a clean, stable regularizing prior
  if (Ncol > 5) {
    for (c in 6:Ncol) {
      betas[c] ~ normal(0, 2.5); // Safe, smooth, handles categoricals, protects Stan engine
    }
  }

  // 3. TiTE Likelihood
  vector[Ntotal] eta = X * betas;
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      target += log_inv_logit(eta[i]);
    } else {
      target += w[i] * log1m_inv_logit(eta[i]);
    }
  }
}

generated quantities {
  // Optional block to cleanly capture log-likelihood for WAIC/LOO diagnostic checks
  vector[Ntotal] log_lik;
  vector[Ntotal] eta_gq = X * betas;
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      log_lik[i] = log_inv_logit(eta_gq[i]);
    } else {
      log_lik[i] = w[i] * log1m_inv_logit(eta_gq[i]);
    }
  }
}