// Name: BioTiTE_V3.stan
// Auth: u.niazi@soton.ac.uk
// Date: 09/09/2026
// Desc: TiTE-CRM Model with dynamic skeleton dose-toxicity priors, 
//       maintaining dose monotonicity and hierarchical tau-shrinkage 
//       for genomic biomarker screening across sequential trial cohorts.

data {
  int<lower=1> Ntotal;                   // Total patients observed
  int<lower=1> Ncol;                     // Total columns in design matrix (5 dose + P covariates)
  matrix[Ntotal, Ncol] X;                // Design matrix
  array[Ntotal] int<lower=0, upper=1> y; // Binary DLT outcomes
  vector[Ntotal] w;                      // TiTE partial follow-up weights
  
  // USER-FRIENDLY SKELETON INPUTS (Full 5-Dose Vectors)
  vector[5] prior_means;                 // Skeleton logit means for Doses 1 to 5
  real<lower=0> prior_sd_alpha1;          // Prior SD for Dose 1 baseline
  vector<lower=0>[4] prior_sds_delta;     // Prior SDs for dose increments (Doses 2 to 5)
}

transformed data {
  int P = Ncol - 5;                      // Number of genomic covariates
  
  // Automatically calculate logit increments inside Stan to protect User Experience
  vector[4] prior_means_delta;
  for (d in 1:4) {
    prior_means_delta[d] = prior_means[d + 1] - prior_means[d];
  }
}

parameters {
  real alpha_1;                          // Baseline logit-risk for Dose 1
  vector<lower=0>[4] delta;              // Strictly positive dose increments (Doses 2 to 5)
  
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

  // 2. Non-centered scaling
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
  // Dynamic clinical skeleton priors using internally calculated deltas
  alpha_1 ~ normal(prior_means[1], prior_sd_alpha1);
  delta ~ normal(prior_means_delta, prior_sds_delta);

  // Global genomic shrinkage priors
  // Loosened hierarchical tau prior (Marginal SD ≈ 2.83)
  tau ~ exponential(0.5);
  beta_cov_raw ~ normal(0, 1.0);

  // Likelihood
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
