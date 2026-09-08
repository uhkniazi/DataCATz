# File: 03_BioTiTE_Skeleton_Sensitivity.R
# Auth: u.niazi@soton.ac.uk
# Date: 08/09/2026
# Desc: Testing various priors i.e. relaxed to strong for doses

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

# ------------------------------------------------------------------------------
# 1. PARAMETER SETUP & SCENARIO 2 DATA GENERATION
# ------------------------------------------------------------------------------
n       <- 40
n_doses <- 5
n_genes <- 5

# Target Skeleton & Logit Intercepts
skeleton   <- c(0.05, 0.10, 0.20, 0.35, 0.50)
alpha_true <- qlogis(skeleton)

# Ground Truth Gene Effects (2 Causal, 3 Null)
beta_true <- c(
  Gene1 = 1.0,  # Strong causal (Confounding target)
  Gene2 = 0.5,  # Moderate causal
  Gene3 = 0.0,  # Null
  Gene4 = 0.0,  # Null
  Gene5 = 0.0   # Null
)

# Phase I Dose Allocation (Skewed to lower/mid doses)
dose_probs <- c(0.30, 0.30, 0.20, 0.10, 0.10)
dose       <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)

# --- SCENARIO 2: Dose-Dependent Gene1 Prevalence ---
# Gene1 probability increases with assigned dose cohort
g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)

G <- matrix(rnorm(n * n_genes), nrow = n, ncol = n_genes)
G[, 1] <- rbinom(n, size = 1, prob = g1_prob[dose]) # Confounded binary trait
colnames(G) <- paste0("Gene", 1:n_genes)

# Linear Predictor & DLT Outcome Generation
eta <- alpha_true[dose] + as.vector(G %*% beta_true)
p   <- plogis(eta)
dlt <- rbinom(n, size = 1, prob = p)

df_sim <- data.frame(
  patient_id = 1:n,
  dose       = dose,
  G,
  p_true     = round(p, 3),
  dlt        = dlt
)

cat("=== SCENARIO 2 DATA SUMMARY ===\n")
cat("Dose Distribution & Observed DLTs:\n")
print(table(Dose = df_sim$dose, DLT = df_sim$dlt))

cat("\nGene 1 Prevalence by Assigned Dose Level:\n")
print(round(tapply(df_sim$Gene1, df_sim$dose, mean), 2))

# ------------------------------------------------------------------------------
# 2. DESIGN MATRIX & STAN DATA ASSEMBLY
# ------------------------------------------------------------------------------
X_design <- model.matrix(~ 0 + factor(dose, levels = 1:5) + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, data = df_sim)

stan_data <- list(
  Ntotal      = nrow(X_design),
  Ncol        = ncol(X_design),
  X           = X_design,
  y           = as.array(df_sim$dlt),
  w           = rep(1.0, nrow(df_sim)),
  prior_means = alpha_true
)

# ==============================================================================
# BIOTITE PRIOR SENSITIVITY TEST: LOOSE vs MODERATE vs STRONG ANCHORING
# ==============================================================================

# 1. Parameterized Stan Model Code
stan_code_sensitivity <- "
data {
  int<lower=1> Ntotal;
  int<lower=1> Ncol;
  matrix[Ntotal, Ncol] X;
  array[Ntotal] int<lower=0, upper=1> y;
  vector[Ntotal] w;
  vector[5] prior_means;
  real<lower=0> sd_alpha;  // Prior SD for baseline Dose 1
  real<lower=0> sd_delta;  // Prior SD for dose increments
}
parameters {
  real alpha_1;
  vector<lower=0>[4] delta;
  vector[Ncol - 5] beta_cov;
}
transformed parameters {
  vector[5] alpha;
  vector[Ntotal] eta;

  alpha[1] = alpha_1;
  alpha[2] = alpha[1] + delta[1];
  alpha[3] = alpha[2] + delta[2];
  alpha[4] = alpha[3] + delta[3];
  alpha[5] = alpha[4] + delta[4];

  for (i in 1:Ntotal) {
    real dose_intercept = 0;
    for (d in 1:5) {
      if (X[i, d] == 1) dose_intercept = alpha[d];
    }
    if (Ncol > 5) {
      vector[Ncol - 5] patient_covariates = to_vector(X[i, 6:Ncol]);
      eta[i] = dose_intercept + dot_product(patient_covariates, beta_cov);
    } else {
      eta[i] = dose_intercept;
    }
  }
}
model {
  alpha_1 ~ normal(prior_means[1], sd_alpha);
  
  delta[1] ~ normal(prior_means[2] - prior_means[1], sd_delta);
  delta[2] ~ normal(prior_means[3] - prior_means[2], sd_delta);
  delta[3] ~ normal(prior_means[4] - prior_means[3], sd_delta);
  delta[4] ~ normal(prior_means[5] - prior_means[4], sd_delta);
  
  beta_cov ~ normal(0, 2.5);

  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      target += log_inv_logit(eta[i]);
    } else {
      target += w[i] * log1m_inv_logit(eta[i]);
    }
  }
}
"

compiled_sens <- rstan::stan_model(model_code = stan_code_sensitivity)

# 2. Fit Across the 3 Prior Configurations
run_biotite_tier <- function(sd_a, sd_d) {
  s_data <- stan_data
  s_data$sd_alpha <- sd_a
  s_data$sd_delta <- sd_d
  
  fit <- rstan::sampling(compiled_sens, data = s_data, iter = 2000, warmup = 1000, chains = 4, refresh = 0)
  ext <- rstan::extract(fit)
  
  genes <- t(apply(ext$beta_cov, 2, function(x) c(Mean = mean(x), SD = sd(x))))
  doses <- apply(plogis(ext$alpha), 2, mean)
  return(list(genes = genes, doses = doses))
}

m1_loose    <- run_biotite_tier(1.50, 1.00)
m2_moderate <- run_biotite_tier(1.00, 0.50)
m3_strong   <- run_biotite_tier(0.50, 0.25)

# 3. Print Comparative Output Matrix
comp_genes <- data.frame(
  True_Effect   = beta_true,
  Loose_Mean    = round(m1_loose$genes[, "Mean"], 3),
  Loose_SD      = round(m1_loose$genes[, "SD"], 3),
  Mod_Mean      = round(m2_moderate$genes[, "Mean"], 3),
  Mod_SD        = round(m2_moderate$genes[, "SD"], 3),
  Strong_Mean   = round(m3_strong$genes[, "Mean"], 3),
  Strong_SD     = round(m3_strong$genes[, "SD"], 3)
)

comp_doses <- data.frame(
  Skeleton      = skeleton,
  Loose_Prob    = round(m1_loose$doses, 3),
  Mod_Prob      = round(m2_moderate$doses, 3),
  Strong_Prob   = round(m3_strong$doses, 3)
)

cat("\n=== SENSITIVITY RESULTS: DOSE CURVE SHIFT ===\n")
print(comp_doses)

cat("\n=== SENSITIVITY RESULTS: GENE COEFFICIENT RECOVERY ===\n")
print(comp_genes)
