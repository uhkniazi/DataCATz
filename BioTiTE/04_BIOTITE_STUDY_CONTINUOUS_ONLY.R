# ==============================================================================
# File: 04_BIOTITE_STUDY_CONTINUOUS_ONLY.R
# Auth: u.niazi@soton.ac.uk
# Date: 09/09/2026
# Desc: High-dimensional ($P=100$) continuous feature debug script. Evaluates 
#       biomarker prioritization and dose-unconfounding using BioTiTE_V3.stan
#       when all predictors are continuous Gaussian variables.
# ==============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

# ==============================================================================
# 1. HELPER: RANK-BASED ROC-AUC FUNCTION (PURE BASE R)
# ==============================================================================
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1)
  n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  auc <- (sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}

# ==============================================================================
# 2. DATA GENERATION: ALL CONTINUOUS FEATURES (P = 100)
# ==============================================================================
generate_continuous_scenario_data <- function(n = 40, 
                                              n_doses = 5, 
                                              n_features = 100, 
                                              skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                              dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  alpha_true <- qlogis(skeleton)
  dose_raw   <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
  dose_fac   <- factor(dose_raw, levels = 1:n_doses)
  
  # --- Ground Truth Coefficients ---
  beta_true <- rep(0.0, n_features)
  names(beta_true) <- paste0("X_", 1:n_features)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) # 5 Causal continuous signals
  
  # --- Continuous Covariate Matrix ---
  X_mat <- matrix(0, nrow = n, ncol = n_features)
  colnames(X_mat) <- paste0("X_", 1:n_features)
  
  dose_scaled <- as.numeric(scale(dose_raw))
  
  # Feature 1 (Causal): Continuous & Strongly Dose-Confounded
  X_mat[, 1] <- rnorm(n, mean = 1.0 * dose_scaled, sd = 1.0)
  
  # Feature 2 (Causal): Continuous & Moderately Dose-Confounded
  X_mat[, 2] <- rnorm(n, mean = 0.6 * dose_scaled, sd = 1.0)
  
  # Features 3-5 (Causal): Continuous & Unconfounded
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  X_mat[, 4] <- rnorm(n, mean = 0, sd = 1.0)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # Features 6-100 (Null Background Noise): Standard Gaussian
  for (j in 6:n_features) {
    X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
  }
  
  # Linear Predictor & DLT Generation
  eta <- alpha_true[dose_raw] + as.vector(X_mat %*% beta_true)
  p   <- plogis(eta)
  dlt <- rbinom(n, size = 1, prob = p)
  
  df_sim <- data.frame(patient_id = 1:n, dose = dose_fac, dlt = dlt, X_mat)
  X_design <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:n_features))])
  
  return(list(
    df_sim     = df_sim,
    X_design   = X_design,
    beta_true  = beta_true,
    alpha_true = alpha_true,
    skeleton   = skeleton
  ))
}

# ==============================================================================
# 3. RUN DEBUG TRIAL & ASSEMBLE MODEL INPUTS
# ==============================================================================
sim_data  <- generate_continuous_scenario_data(n = 40, n_features = 100)
df_sim    <- sim_data$df_sim
X_design  <- sim_data$X_design
beta_true <- sim_data$beta_true

cat("=== CONTINUOUS DATA SUMMARY ===\n")
cat("Cohort Size:", nrow(df_sim), "| Total Features:", length(beta_true), "\n")
cat("Dose Allocations:\n")
print(table(df_sim$dose))

stan_data_std <- list(
  Ntotal = nrow(X_design),
  Ncol   = ncol(X_design),
  X      = X_design,
  y      = as.array(df_sim$dlt)
)

stan_data_biotite <- list(
  Ntotal          = nrow(X_design),
  Ncol            = ncol(X_design),
  X               = X_design,
  y               = as.array(df_sim$dlt),
  w               = rep(1.0, nrow(df_sim)),
  prior_means     = qlogis(sim_data$skeleton),
  prior_sd_alpha1 = 0.50,
  prior_sds_delta = rep(0.25, 4)
)

# ==============================================================================
# 4. MODEL FITTING (TIER 0, TIER 1, TIER 2)
# ==============================================================================
# --- Tier 0: GLM ---
fit_glm <- suppressWarnings(
  glm(dlt ~ ., data = df_sim[, -1], family = binomial(link = "logit"))
)
glm_coef_mat <- summary(fit_glm)$coefficients

glm_beta <- rep(0, 100)
names(glm_beta) <- paste0("X_", 1:100)
matched_names <- intersect(names(glm_beta), rownames(glm_coef_mat))
glm_beta[matched_names] <- glm_coef_mat[matched_names, "Estimate"]

# --- Tier 1: Standard Stan ---
compiled_std <- rstan::stan_model(file = "binomialRegression.stan")
fit_std <- rstan::sampling(compiled_std, data = stan_data_std, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_std  <- rstan::extract(fit_std)
std_beta <- apply(ext_std$betas[, 6:105], 2, mean)
names(std_beta) <- paste0("X_", 1:100)

# --- Tier 2: BioTiTE V3 ---
compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")
fit_biotite <- rstan::sampling(compiled_biotite, data = stan_data_biotite, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_bio  <- rstan::extract(fit_biotite)
bio_beta <- apply(ext_bio$beta_cov, 2, mean)
names(bio_beta) <- paste0("X_", 1:100)

# ==============================================================================
# 5. PRIORITIZATION METRICS (RANKS & ROC-AUC)
# ==============================================================================
causal_mask <- ifelse(beta_true > 0, 1, 0)

glm_ranks <- rank(abs(glm_beta))
std_ranks <- rank(abs(std_beta))
bio_ranks <- rank(abs(bio_beta))

rank_summary <- data.frame(
  True_Effect   = beta_true[1:5],
  Is_Confounded = c("Yes (Strong Cont)", "Yes (Mod Cont)", "No", "No", "No"),
  GLM_Rank      = glm_ranks[1:5],
  StdStan_Rank  = std_ranks[1:5],
  BioTiTE_Rank  = bio_ranks[1:5]
)

auc_summary <- data.frame(
  Model = c("GLM (Tier 0)", "Standard Stan (Tier 1)", "BioTiTE V3 (Tier 2)"),
  Biomarker_AUC = c(
    calc_feature_auc(abs(glm_beta), causal_mask),
    calc_feature_auc(abs(std_beta), causal_mask),
    calc_feature_auc(abs(bio_beta), causal_mask)
  )
)

cat("\n======================================================================\n")
cat("TOP CAUSAL FEATURE RANKINGS (All Continuous Features; Higher is better)\n")
cat("======================================================================\n")
print(rank_summary)

cat("\n======================================================================\n")
cat("GLOBAL FEATURE PRIORITIZATION ROC-AUC (Causal vs 95 Null Features)\n")
cat("======================================================================\n")
print(auc_summary)

# ==============================================================================
# 6. EXTRACT SKELETON POSTERIORS FOR SEQUENTIAL LEARNING (TRIAL 1 -> TRIAL 2)
# ==============================================================================
alpha_post_means <- apply(ext_bio$alpha, 2, mean)
alpha_post_sds   <- apply(ext_bio$alpha, 2, sd)
delta_post_sds   <- apply(ext_bio$delta, 2, sd)

cat("\n======================================================================\n")
cat("EXTRACTED DOSE POSTERIORS FOR TRIAL 2 PRIOR INJECTION\n")
cat("======================================================================\n")
post_prior_df <- data.frame(
  Dose_Level           = paste0("Dose_", 1:5),
  Posterior_Mean_Alpha = round(alpha_post_means, 3),
  Posterior_SD_Alpha   = round(alpha_post_sds, 3),
  Fitted_Skeleton_Prob = round(plogis(alpha_post_means), 3)
)
print(post_prior_df)