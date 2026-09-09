// Name: BioTiTE_V2.stan
// Auth: u.niazi@soton.ac.uk
// Date: 08/09/2026
// Desc: TiTE-CRM Model with skeleton dose toxicity priors, maintaining
//       monotonic behaviour of toxicity prediction and using covariates
//       with a Hierarchical adaptive prior 
//       to screen for features associated with toxicity 

data {
  int<lower=1> Ntotal;                   // Total patients observed
  int<lower=1> Ncol;                     // Total columns in model matrix (5 dose + P covariates)
  matrix[Ntotal, Ncol] X;                // Design matrix (cols 1:5 = dose indicators, 6:Ncol = covariates)
  array[Ntotal] int<lower=0, upper=1> y; // Binary DLT outcomes
  vector[Ntotal] w;                      // TiTE partial follow-up weights
  vector[5] prior_means;                 // Logit-scale target skeleton
}

transformed data {
  int P = Ncol - 5;                      // Total number of genomic covariates
}

parameters {
  real alpha_1;                          // Baseline logit-risk for Dose 1
  vector<lower=0>[4] delta;              // Positive dose increments (Doses 2 to 5)
  
  // Non-Centered Genomic Shrinkage Parameters
  vector[P] beta_cov_raw;                // Unscaled covariate effects ~ N(0,1)
  real<lower=0> tau;                     // Global scale parameter
}

transformed parameters {
  vector[5] alpha;
  vector[P] beta_cov;
  vector[Ntotal] eta;

  // 1. Construct strictly monotonic dose-response intercepts
  alpha[1] = alpha_1;
  alpha[2] = alpha[1] + delta[1];
  alpha[3] = alpha[2] + delta[2];
  alpha[4] = alpha[3] + delta[3];
  alpha[5] = alpha[4] + delta[4];

  // 2. Non-centered scaling: eliminates Neal's funnel as tau -> 0
  beta_cov = beta_cov_raw * tau;

  // 3. Assemble linear predictor
  for (i in 1:Ntotal) {
    real dose_intercept = 0;
    for (d in 1:5) {
      if (X[i, d] == 1) dose_intercept = alpha[d];
    }
    
    if (P > 0) {
      vector[P] patient_covariates = to_vector(X[i, 6:Ncol]);
      eta[i] = dose_intercept + dot_product(patient_covariates, beta_cov);
    } else {
      eta[i] = dose_intercept;
    }
  }
}

model {
  // ===========================================================================
  // PRIORS: STRONG SKELETON ANCHOR + HIERARCHICAL GENOMIC SHRINKAGE
  // ===========================================================================
  
  // Strong clinical skeleton anchor (Locks down dose geometry)
  alpha_1 ~ normal(prior_means[1], 0.50);
  delta[1] ~ normal(prior_means[2] - prior_means[1], 0.25);
  delta[2] ~ normal(prior_means[3] - prior_means[2], 0.25);
  delta[3] ~ normal(prior_means[4] - prior_means[3], 0.25);
  delta[4] ~ normal(prior_means[5] - prior_means[4], 0.25);

  // Global genomic shrinkage priors
  tau ~ exponential(1.0);
  beta_cov_raw ~ normal(0, 1.0);

  // ===========================================================================
  // LIKELIHOOD (TiTE Partial Likelihood)
  // ===========================================================================
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      target += log_inv_logit(eta[i]);
    } else {
      target += w[i] * log1m_inv_logit(eta[i]);
    }
  }
}

generated quantities {
  vector[Ntotal] log_lik;
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      log_lik[i] = log_inv_logit(eta[i]);
    } else {
      log_lik[i] = w[i] * log1m_inv_logit(eta[i]);
    }
  }
}
