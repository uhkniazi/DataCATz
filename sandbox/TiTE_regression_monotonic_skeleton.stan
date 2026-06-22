data {
  int<lower=1> Ntotal;                  // Total number of patients observed so far
  int<lower=1> Ncol;                    // Total columns in design matrix (5 dose indicators + covariates)
  matrix[Ntotal, Ncol] X;               // Design matrix (first 5 columns are dose indicators)
  array[Ntotal] int<lower=0, upper=1> y;// Binary DLT outcome (1 = DLT, 0 = No DLT)
  vector[Ntotal] w;                     // TiTE partial follow-up weights (range 0.0 to 1.0)
  vector[5] prior_means;                // Logit-scale skeleton anchoring points
}

parameters {
  real alpha_1;                         // Baseline logit-risk for Dose 1
  vector<lower=0>[4] delta;             // Strictly positive increments for Doses 2 through 5
  vector[Ncol - 5] beta_cov;            // Regression coefficients for patient covariates (Sex, BMI)
}

transformed parameters {
  vector[5] alpha;
  vector[Ntotal] eta;

  // 1. Construct the strictly increasing dose-response curve via sequential addition
  alpha[1] = alpha_1;
  alpha[2] = alpha[1] + delta[1];
  alpha[3] = alpha[2] + delta[2];
  alpha[4] = alpha[3] + delta[3];
  alpha[5] = alpha[4] + delta[4];

  // 2. Compute individual linear predictors
  for (i in 1:Ntotal) {
    real dose_intercept = 0;
    
    // Scan the first 5 columns of the design matrix to find the assigned dose level
    for (d in 1:5) {
      if (X[i, d] == 1) {
        dose_intercept = alpha[d];
      }
    }
    
    // Combine the dose intercept with patient-specific covariate shifts if present
    if (Ncol > 5) {
      vector[Ncol - 5] patient_covariates = to_vector(X[i, 6:Ncol]);
      eta[i] = dose_intercept + dot_product(patient_covariates, beta_cov);
    } else {
      eta[i] = dose_intercept;
    }
  }
}

model {
  // ===========================================================================
  // PRIORS (Relaxed to sd = 1.0/1.5 to prioritize data-driven curve bending)
  // ===========================================================================
  
  // Anchor the baseline Dose 1 risk near the skeleton mean, but allow movement
  alpha_1 ~ normal(prior_means[1], 1.5);
  
  // Relaxed independent priors on the gaps (centered on skeleton increments)
  delta[1] ~ normal(prior_means[2] - prior_means[1], 1.0); 
  delta[2] ~ normal(prior_means[3] - prior_means[2], 1.0);
  delta[3] ~ normal(prior_means[4] - prior_means[3], 1.0);
  delta[4] ~ normal(prior_means[5] - prior_means[4], 1.0);
  
  // Weakly informative prior for patient traits (Sex, BMI, etc.)
  beta_cov ~ normal(0, 2.5);

  // ===========================================================================
  // LIKELIHOOD (Time-to-Event Partial Likelihood)
  // ===========================================================================
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      // If a DLT occurs, use the full log-probability of the event
      target += log_inv_logit(eta[i]);
    } else {
      // If no DLT has occurred yet, scale survival probability by follow-up time
      target += w[i] * log1m_inv_logit(eta[i]);
    }
  }
}

generated quantities {
  vector[Ntotal] log_lik;
  
  // Retain point-by-point log-likelihoods for WAIC/LOO diagnostic metrics
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      log_lik[i] = log_inv_logit(eta[i]);
    } else {
      log_lik[i] = w[i] * log1m_inv_logit(eta[i]);
    }
  }
}
