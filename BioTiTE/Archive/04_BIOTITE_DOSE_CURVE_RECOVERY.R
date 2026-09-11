# ==============================================================================
# File: 04_BIOTITE_DOSE_CURVE_RECOVERY.R
# Auth: u.niazi@soton.ac.uk
# Date: 09/09/2026
# Desc: Evaluates dose-toxicity curve recovery (MSE & Monotonicity) under 
#       prior skeleton misspecification and P=100 covariate adjustment.
# ==============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

# ==============================================================================
# 1. DATA GENERATION WITH SKELETON MISSPECIFICATION
# ==============================================================================
generate_curve_recovery_data <- function(n = 40, 
                                         n_doses = 5, 
                                         n_features = 100, 
                                         true_skeleton  = c(0.08, 0.15, 0.28, 0.45, 0.62),
                                         dose_probs     = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  alpha_true <- qlogis(true_skeleton)
  dose_raw   <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
  dose_fac   <- factor(dose_raw, levels = 1:n_doses)
  
  # ground truth effects (5 causal, 95 nulls)
  beta_true <- rep(0.0, n_features)
  names(beta_true) <- paste0("X_", 1:n_features)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) 
  
  X_mat <- matrix(0, nrow = n, ncol = n_features)
  colnames(X_mat) <- paste0("X_", 1:n_features)
  
  # Feature 1: Centered Binary {-0.5, +0.5} & Dose-Confounded
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  raw_b1  <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  X_mat[, 1] <- ifelse(raw_b1 == 1, 0.5, -0.5)
  
  # Feature 2: Continuous & Dose-Confounded
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  
  # Features 3-5: Unconfounded signals
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  raw_b4  <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 4] <- ifelse(raw_b4 == 1, 0.5, -0.5)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # Features 6-100: Null background noise
  for (j in 6:n_features) {
    if (j %% 2 == 0) {
      X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
    } else {
      raw_null <- rbinom(n, size = 1, prob = 0.30)
      X_mat[, j] <- ifelse(raw_null == 1, 0.5, -0.5)
    }
  }
  
  # Generative Linear Predictor
  eta <- alpha_true[dose_raw] + as.vector(X_mat %*% beta_true)
  p   <- plogis(eta)
  dlt <- rbinom(n, size = 1, prob = p)
  
  df_sim <- data.frame(patient_id = 1:n, dose = dose_fac, dlt = dlt, X_mat)
  
  # Design matrix WITH covariates (105 columns)
  X_design_cov <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:n_features))])
  
  # Design matrix WITHOUT covariates (5 dose columns only)
  X_design_nocov <- model.matrix(~ 0 + dose, data = df_sim)
  
  return(list(
    df_sim         = df_sim,
    X_design_cov   = X_design_cov,
    X_design_nocov = X_design_nocov,
    true_skeleton  = true_skeleton,
    beta_true      = beta_true
  ))
}

# ==============================================================================
# 2. RUN DATA GENERATION
# ==============================================================================
# Misspecified clinical prior given to BioTiTE model
prior_skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)

sim_data <- generate_curve_recovery_data(n = 40, true_skeleton = c(0.08, 0.15, 0.28, 0.45, 0.62))
df_sim   <- sim_data$df_sim

cat("=== DOSE CURVE RECOVERY EXPERIMENT ===\n")
cat("Ground Truth Curve:  ", sim_data$true_skeleton, "\n")
cat("BioTiTE Prior Input: ", prior_skeleton, "\n")
cat("Dose Allocations:\n")
print(table(df_sim$dose))

# ==============================================================================
# 3. HELPER FUNCTIONS FOR CURVE EVALUATION
# ==============================================================================
calc_curve_mse <- function(est_probs, true_probs) {
  mean((est_probs - true_probs)^2)
}

check_mono_violation <- function(est_probs) {
  any(diff(est_probs) < 0)
}

# ==============================================================================
# 4. MODEL FITTING ACROSS ALL 6 CONFIGURATIONS
# ==============================================================================
compiled_std <- rstan::stan_model(file = "binomialRegression.stan")
compiled_bio <- rstan::stan_model(file = "BioTiTE_V3.stan")

# --- A. GLM (No Covariates) ---
fit_glm_nocov  <- glm(dlt ~ 0 + dose, data = df_sim, family = binomial(link = "logit"))
p_glm_nocov    <- plogis(coef(fit_glm_nocov)[1:5])

# --- B. GLM (With Covariates) ---
fit_glm_cov    <- suppressWarnings(glm(dlt ~ ., data = df_sim[, -1], family = binomial(link = "logit")))
coef_mat_glm   <- summary(fit_glm_cov)$coefficients
p_glm_cov      <- rep(NA, 5)
for(d in 1:5) {
  row_name <- paste0("dose", d)
  if(row_name %in% rownames(coef_mat_glm)) p_glm_cov[d] <- plogis(coef_mat_glm[row_name, "Estimate"])
}

# --- C. Standard Stan (No Covariates) ---
stan_data_std_nocov <- list(Ntotal = nrow(sim_data$X_design_nocov), Ncol = 5, X = sim_data$X_design_nocov, y = as.array(df_sim$dlt))
fit_std_nocov <- rstan::sampling(compiled_std, data = stan_data_std_nocov, iter = 2000, warmup = 1000, chains = 4, refresh = 0)
p_std_nocov   <- apply(plogis(rstan::extract(fit_std_nocov)$betas[, 1:5]), 2, mean)

# --- D. Standard Stan (With Covariates) ---
stan_data_std_cov <- list(Ntotal = nrow(sim_data$X_design_cov), Ncol = 105, X = sim_data$X_design_cov, y = as.array(df_sim$dlt))
fit_std_cov <- rstan::sampling(compiled_std, data = stan_data_std_cov, iter = 2000, warmup = 1000, chains = 4, refresh = 0)
p_std_cov   <- apply(plogis(rstan::extract(fit_std_cov)$betas[, 1:5]), 2, mean)

# --- E. BioTiTE V3 (No Covariates) ---
stan_data_bio_nocov <- list(
  Ntotal = nrow(sim_data$X_design_nocov), Ncol = 5, X = sim_data$X_design_nocov, y = as.array(df_sim$dlt), w = rep(1.0, 40),
  prior_means = qlogis(prior_skeleton), prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
)
fit_bio_nocov <- rstan::sampling(compiled_bio, data = stan_data_bio_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
p_bio_nocov   <- apply(plogis(rstan::extract(fit_bio_nocov)$alpha), 2, mean)

# --- F. BioTiTE V3 (With Covariates) ---
stan_data_bio_cov <- list(
  Ntotal = nrow(sim_data$X_design_cov), Ncol = 105, X = sim_data$X_design_cov, y = as.array(df_sim$dlt), w = rep(1.0, 40),
  prior_means = qlogis(prior_skeleton), prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
)
fit_bio_cov <- rstan::sampling(compiled_bio, data = stan_data_bio_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
p_bio_cov   <- apply(plogis(rstan::extract(fit_bio_cov)$alpha), 2, mean)

# ==============================================================================
# 5. DIAGNOSTIC COMPARISON & CURVE RECOVERY METRICS
# ==============================================================================
results_df <- data.frame(
  Model_Tier          = c("GLM (Tier 0)", "GLM (Tier 0)", 
                          "Standard Stan (Tier 1)", "Standard Stan (Tier 1)", 
                          "BioTiTE V3 (Tier 2)", "BioTiTE V3 (Tier 2)"),
  Covariates_Fitted   = rep(c("No (Omitted)", "Yes (P=100)"), 3),
  Dose_1_Prob         = round(c(p_glm_nocov[1], p_glm_cov[1], p_std_nocov[1], p_std_cov[1], p_bio_nocov[1], p_bio_cov[1]), 3),
  Dose_2_Prob         = round(c(p_glm_nocov[2], p_glm_cov[2], p_std_nocov[2], p_std_cov[2], p_bio_nocov[2], p_bio_cov[2]), 3),
  Dose_3_Prob         = round(c(p_glm_nocov[3], p_glm_cov[3], p_std_nocov[3], p_std_cov[3], p_bio_nocov[3], p_bio_cov[3]), 3),
  Dose_4_Prob         = round(c(p_glm_nocov[4], p_glm_cov[4], p_std_nocov[4], p_std_cov[4], p_bio_nocov[4], p_bio_cov[4]), 3),
  Dose_5_Prob         = round(c(p_glm_nocov[5], p_glm_cov[5], p_std_nocov[5], p_std_cov[5], p_bio_nocov[5], p_bio_cov[5]), 3),
  Curve_MSE           = round(c(
    calc_curve_mse(p_glm_nocov, sim_data$true_skeleton),
    calc_curve_mse(p_glm_cov, sim_data$true_skeleton),
    calc_curve_mse(p_std_nocov, sim_data$true_skeleton),
    calc_curve_mse(p_std_cov, sim_data$true_skeleton),
    calc_curve_mse(p_bio_nocov, sim_data$true_skeleton),
    calc_curve_mse(p_bio_cov, sim_data$true_skeleton)
  ), 4),
  Monotonicity_Viol   = c(
    check_mono_violation(p_glm_nocov),
    check_mono_violation(p_glm_cov),
    check_mono_violation(p_std_nocov),
    check_mono_violation(p_std_cov),
    check_mono_violation(p_bio_nocov),
    check_mono_violation(p_bio_cov)
  )
)

cat("\n======================================================================\n")
cat("DOSE-TOXICITY CURVE RECOVERY DIAGNOSTIC TABLE\n")
cat("Ground Truth Curve: ", sim_data$true_skeleton, "\n")
cat("Prior Skeleton Input: ", prior_skeleton, "\n")
cat("======================================================================\n")
print(results_df)