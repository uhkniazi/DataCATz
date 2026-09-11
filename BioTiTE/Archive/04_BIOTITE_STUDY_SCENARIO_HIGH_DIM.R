# ==============================================================================
# File: 04_BIOTITE_STUDY_SCENARIO_HIGH_DIM.R
# Auth: u.niazi@soton.ac.uk
# Date: 09/09/2026
# Desc: High-dimensional ($P=100$) debug script evaluating biomarker prioritization,
#       dose-confounding adjustment, and posterior extraction using BioTiTE_V3.stan.
# ==============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

# ==============================================================================
# 1. HELPER: RANK-BASED ROC-AUC FUNCTION (PURE BASE R)
# ==============================================================================
calc_feature_auc <- function(scores, true_is_causal) {
  # scores: absolute magnitude of estimated betas |beta_hat|
  # true_is_causal: binary vector (1 = causal feature, 0 = null feature)
  n1 <- sum(true_is_causal == 1)
  n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  auc <- (sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}

# ==============================================================================
# 2. DATA GENERATION: P = 100 BIOMARKERS (5 CAUSAL, 2 DOSE-CONFOUNDED)
# ==============================================================================
generate_highdim_scenario_data <- function(n = 40, 
                                           n_doses = 5, 
                                           n_features = 100, 
                                           skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                           dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  alpha_true <- qlogis(skeleton)
  dose_raw   <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
  dose_fac   <- factor(dose_raw, levels = 1:n_doses)
  
  # Ground Truth Coefficients: 5 Causal signals, 95 Null signals
  beta_true <- rep(0.0, n_features)
  names(beta_true) <- paste0("X_", 1:n_features)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) 
  
  # Construct Covariate Matrix
  X_mat <- matrix(0, nrow = n, ncol = n_features)
  colnames(X_mat) <- paste0("X_", 1:n_features)
  
  # Feature 1 (Causal): Binary & Dose-Confounded
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  X_mat[, 1] <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  
  # Feature 2 (Causal): Continuous & Dose-Confounded
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  
  # Features 3-5 (Causal): Unconfounded signals
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  X_mat[, 4] <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # Features 6-100 (Null Noise): 50% Gaussian, 50% Bernoulli
  for (j in 6:n_features) {
    if (j %% 2 == 0) {
      X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
    } else {
      X_mat[, j] <- rbinom(n, size = 1, prob = 0.30)
    }
  }
  
  # Outcome generation
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
sim_data  <- generate_highdim_scenario_data(n = 40, n_features = 100)
df_sim    <- sim_data$df_sim
X_design  <- sim_data$X_design
beta_true <- sim_data$beta_true

cat("=== DATA SUMMARY ===\n")
cat("Cohort Size:", nrow(df_sim), "| Total Features:", length(beta_true), "\n")
cat("Dose Allocations:\n")
print(table(df_sim$dose))

stan_data_std <- list(
  Ntotal = nrow(X_design),
  Ncol   = ncol(X_design),
  X      = X_design,
  y      = as.array(df_sim$dlt)
)

# Refactored for BioTiTE_V3.stan dynamic hyperparameter block
stan_data_biotite <- list(
  Ntotal          = nrow(X_design),
  Ncol            = ncol(X_design),
  X               = X_design,
  y               = as.array(df_sim$dlt),
  w               = rep(1.0, nrow(df_sim)),
  prior_means     = qlogis(sim_data$skeleton), # Length 5 vector
  prior_sd_alpha1 = 0.50,                      # Elicited baseline SD
  prior_sds_delta = rep(0.25, 4)               # Elicited increment SDs
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
fit_std <- rstan::sampling(compiled_std, data = stan_data_std, iter = 2000, warmup = 1000, chains = 4, refresh = 0)
ext_std  <- rstan::extract(fit_std)
std_beta <- apply(ext_std$betas[, 6:105], 2, mean)
names(std_beta) <- paste0("X_", 1:100)

# --- Tier 2: BioTiTE V3 ---
compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")
fit_biotite <- rstan::sampling(compiled_biotite, data = stan_data_biotite, iter = 2000, warmup = 1000, chains = 4, refresh = 0)
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
  Is_Confounded = c("Yes (Binary)", "Yes (Cont)", "No", "No", "No"),
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
cat("TOP CAUSAL FEATURE RANKINGS (Out of 100 features; Higher is better)\n")
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

cat("\nExtracted Delta SDs (Increments 1->2, 2->3, 3->4, 4->5):\n")
print(round(delta_post_sds, 3))
