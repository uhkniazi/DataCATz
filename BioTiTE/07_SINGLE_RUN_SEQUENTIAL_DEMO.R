# ==============================================================================
# File: 07_SINGLE_RUN_SEQUENTIAL_DEMO.R
# Auth: u.niazi@soton.ac.uk
# Date: 10/09/2026
# Desc: Single deterministic trial run demonstrating sequential posterior transfer 
#       from Trial 1 to Trial 2 using BioTiTE V3 and Standard Stan GLM.
# ==============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
set.seed(42)

# ------------------------------------------------------------------------------
# A. HELPERS & DATA GENERATOR
# ------------------------------------------------------------------------------
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1); n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  return((sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}

calc_topk_recall <- function(beta_scores, true_causal_idx = 1:5, k = 10) {
  top_k_indices <- order(abs(beta_scores), decreasing = TRUE)[1:k]
  return(sum(top_k_indices %in% true_causal_idx) / length(true_causal_idx))
}

generate_trial_data <- function(n = 40, true_skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50)) {
  alpha_true <- qlogis(true_skeleton)
  
  repeat {
    dose_raw <- sample(1:5, size = n, replace = TRUE, prob = c(0.30, 0.30, 0.20, 0.10, 0.10))
    if (length(unique(dose_raw)) == 5) break
  }
  dose_fac <- factor(dose_raw, levels = 1:5)
  
  beta_true <- rep(0.0, 100)
  names(beta_true) <- paste0("X_", 1:100)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) 
  
  X_mat <- matrix(0, nrow = n, ncol = 100)
  colnames(X_mat) <- paste0("X_", 1:100)
  
  # Features 1-5: True Causal
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  raw_b1  <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  X_mat[, 1] <- ifelse(raw_b1 == 1, 0.5, -0.5)
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  raw_b4  <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 4] <- ifelse(raw_b4 == 1, 0.5, -0.5)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # Features 6-100: Noise
  for (j in 6:100) {
    if (j %% 2 == 0) X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
    else X_mat[, j] <- ifelse(rbinom(n, size = 1, prob = 0.30) == 1, 0.5, -0.5)
  }
  
  p   <- plogis(alpha_true[dose_raw] + as.vector(X_mat %*% beta_true))
  dlt <- rbinom(n, size = 1, prob = p)
  
  df_sim   <- data.frame(patient_id = 1:n, dose = dose_fac, dlt = dlt, X_mat)
  X_design <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:100))])
  
  return(list(df_sim = df_sim, X_design = X_design, true_skeleton = true_skeleton))
}

# ------------------------------------------------------------------------------
# B. PRE-COMPILE STAN MODELS
# ------------------------------------------------------------------------------
cat("Compiling Stan models...\n")
compiled_std <- rstan::stan_model(file = "binomialRegression.stan")
compiled_bio <- rstan::stan_model(file = "BioTiTE_V3.stan")

# Define Skeletons
true_skeleton         <- c(0.05, 0.10, 0.20, 0.35, 0.50)
misspecified_skeleton <- c(0.02, 0.05, 0.10, 0.20, 0.35) # Wrong starting clinical guess

causal_mask <- c(rep(1, 5), rep(0, 95))

# Generate Data
t1_data <- generate_trial_data(n = 40, true_skeleton = true_skeleton)
t2_data <- generate_trial_data(n = 40, true_skeleton = true_skeleton)

# ------------------------------------------------------------------------------
# STEP 1: FIT TRIAL 1 MODELS
# ------------------------------------------------------------------------------
cat("\n=== STEP 1: FITTING TRIAL 1 (N = 40) ===\n")

# A. Standard Stan Trial 1
std_t1_data <- list(Ntotal = 40, Ncol = 105, X = t1_data$X_design, y = as.array(t1_data$df_sim$dlt))
fit_std_t1  <- rstan::sampling(compiled_std, data = std_t1_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_std_t1  <- rstan::extract(fit_std_t1)
p_std_t1    <- apply(plogis(ext_std_t1$betas[, 1:5]), 2, mean)

# B. BioTiTE V3 Trial 1 (Using Intentionally Misspecified Skeleton)
bio_t1_data <- list(
  Ntotal = 40, Ncol = 105, X = t1_data$X_design, y = as.array(t1_data$df_sim$dlt),
  w = rep(1.0, 40), prior_means = qlogis(misspecified_skeleton), 
  prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
)
fit_bio_t1 <- rstan::sampling(compiled_bio, data = bio_t1_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_bio_t1 <- rstan::extract(fit_bio_t1)
p_bio_t1   <- apply(plogis(ext_bio_t1$alpha), 2, mean)

# ------------------------------------------------------------------------------
# STEP 2: EXTRACT POSTERIOR FROM TRIAL 1 & CONSTRUCT TRIAL 2 PRIORS
# ------------------------------------------------------------------------------
cat("\n=== STEP 2: EXTRACTING TRIAL 1 POSTERIOR & SCALING PRIORS ===\n")

t1_alpha_samples <- ext_bio_t1$alpha                              # Draws x 5 matrix
t1_delta_samples <- t1_alpha_samples[, 2:5] - t1_alpha_samples[, 1:4] # Joint delta draws

t1_alpha1_mean <- mean(t1_alpha_samples[, 1])
t1_alpha1_sd   <- sd(t1_alpha_samples[, 1])

t1_delta_means <- colMeans(t1_delta_samples)
t1_delta_sds   <- apply(t1_delta_samples, 2, sd)

# Power Prior Discount Factor (alpha_power = 0.5 -> double the variance)
alpha_power <- 0.5
t2_prior_sd_alpha1   <- t1_alpha1_sd / sqrt(alpha_power)
t2_prior_sds_delta   <- t1_delta_sds / sqrt(alpha_power)
t2_prior_means_alpha <- c(t1_alpha1_mean, t1_alpha1_mean + cumsum(t1_delta_means))

cat("Learned Prior Means for Trial 2 (Logit Scale): ", round(t2_prior_means_alpha, 3), "\n")
cat("Learned Prior Means for Trial 2 (Prob Scale):  ", round(plogis(t2_prior_means_alpha), 3), "\n")
cat("Power-Scaled Prior SDs for Deltas:             ", round(t2_prior_sds_delta, 3), "\n")

# ------------------------------------------------------------------------------
# STEP 3: FIT TRIAL 2 MODELS
# ------------------------------------------------------------------------------
cat("\n=== STEP 3: FITTING TRIAL 2 (N = 40) ===\n")

# A. Standard Stan Trial 2 (Independent fit)
std_t2_data <- list(Ntotal = 40, Ncol = 105, X = t2_data$X_design, y = as.array(t2_data$df_sim$dlt))
fit_std_t2  <- rstan::sampling(compiled_std, data = std_t2_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_std_t2  <- rstan::extract(fit_std_t2)
p_std_t2    <- apply(plogis(ext_std_t2$betas[, 1:5]), 2, mean)

# B. BioTiTE V3 Trial 2 (FRESH PRIOR - Ignores Trial 1)
bio_t2_fresh_data <- list(
  Ntotal = 40, Ncol = 105, X = t2_data$X_design, y = as.array(t2_data$df_sim$dlt),
  w = rep(1.0, 40), prior_means = qlogis(misspecified_skeleton), 
  prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
)
fit_bio_t2_fresh <- rstan::sampling(compiled_bio, data = bio_t2_fresh_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_bio_t2_fresh <- rstan::extract(fit_bio_t2_fresh)
p_bio_t2_fresh   <- apply(plogis(ext_bio_t2_fresh$alpha), 2, mean)

# C. BioTiTE V3 Trial 2 (SEQUENTIAL PRIOR - Inherits Trial 1 Posterior)
bio_t2_seq_data <- list(
  Ntotal = 40, Ncol = 105, X = t2_data$X_design, y = as.array(t2_data$df_sim$dlt),
  w = rep(1.0, 40), prior_means = t2_prior_means_alpha, 
  prior_sd_alpha1 = t2_prior_sd_alpha1, prior_sds_delta = t2_prior_sds_delta
)
fit_bio_t2_seq <- rstan::sampling(compiled_bio, data = bio_t2_seq_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
ext_bio_t2_seq <- rstan::extract(fit_bio_t2_seq)
p_bio_t2_seq   <- apply(plogis(ext_bio_t2_seq$alpha), 2, mean)

# ------------------------------------------------------------------------------
# STEP 4: COMPUTE METRICS & DISPLAY RESULTS
# ------------------------------------------------------------------------------
rmse <- function(est, true) sqrt(mean((est - true)^2))

bio_t1_beta     <- apply(ext_bio_t1$beta_cov, 2, mean)
bio_t2_f_beta   <- apply(ext_bio_t2_fresh$beta_cov, 2, mean)
bio_t2_s_beta   <- apply(ext_bio_t2_seq$beta_cov, 2, mean)

curve_table <- data.frame(
  Dose_Level         = paste0("Dose_", 1:5),
  True_Skeleton      = true_skeleton,
  Misspecified_Input = misspecified_skeleton,
  StdStan_T1         = round(p_std_t1, 3),
  BioTiTE_T1         = round(p_bio_t1, 3),
  StdStan_T2         = round(p_std_t2, 3),
  BioTiTE_T2_Fresh   = round(p_bio_t2_fresh, 3),
  BioTiTE_T2_Seq     = round(p_bio_t2_seq, 3)
)

rmse_table <- data.frame(
  Model_Fit                 = c("Std Stan T1", "BioTiTE T1 (Misspecified Input)", 
                                "Std Stan T2", "BioTiTE T2 (Fresh Input)", 
                                "BioTiTE T2 (Sequential Update)"),
  Dose_Curve_RMSE           = round(c(rmse(p_std_t1, true_skeleton), 
                                      rmse(p_bio_t1, true_skeleton), 
                                      rmse(p_std_t2, true_skeleton), 
                                      rmse(p_bio_t2_fresh, true_skeleton), 
                                      rmse(p_bio_t2_seq, true_skeleton)), 4),
  Biomarker_AUC             = round(c(calc_feature_auc(abs(apply(ext_std_t1$betas[, 6:105], 2, mean)), causal_mask),
                                      calc_feature_auc(abs(bio_t1_beta), causal_mask),
                                      calc_feature_auc(abs(apply(ext_std_t2$betas[, 6:105], 2, mean)), causal_mask),
                                      calc_feature_auc(abs(bio_t2_f_beta), causal_mask),
                                      calc_feature_auc(abs(bio_t2_s_beta), causal_mask)), 4),
  Biomarker_Top10_Recall    = round(c(calc_topk_recall(abs(apply(ext_std_t1$betas[, 6:105], 2, mean)), 1:5, k = 10),
                                      calc_topk_recall(abs(bio_t1_beta), 1:5, k = 10),
                                      calc_topk_recall(abs(apply(ext_std_t2$betas[, 6:105], 2, mean)), 1:5, k = 10),
                                      calc_topk_recall(abs(bio_t2_f_beta), 1:5, k = 10),
                                      calc_topk_recall(abs(bio_t2_s_beta), 1:5, k = 10)), 4)
)

cat("\n======================================================================\n")
cat("ESTIMATED DOSE-TOXICITY PROBABILITIES (SINGLE TRIAL DEMO)\n")
cat("======================================================================\n")
print(curve_table)

cat("\n======================================================================\n")
cat("MODEL PERFORMANCE COMPARISON TABLE\n")
cat("======================================================================\n")
print(rmse_table)

######### added example of sequential borrowing sweep
# Test borrowing strength discounting: alpha_power from 0.05 (weak) to 1.0 (strong)
alpha_powers <- c(0.05, 0.10, 0.25, 0.50, 1.00)
sweep_results <- data.frame()

for (ap in alpha_powers) {
  # Apply power-prior scaling
  t2_sd_a1 <- t1_alpha1_sd / sqrt(ap)
  t2_sds_d <- t1_delta_sds / sqrt(ap)
  
  data_seq <- list(
    Ntotal = 40, Ncol = 105, X = t2_data$X_design, y = as.array(t2_data$df_sim$dlt),
    w = rep(1.0, 40), prior_means = t2_prior_means_alpha, 
    prior_sd_alpha1 = t2_sd_a1, prior_sds_delta = t2_sds_d
  )
  
  fit_seq <- rstan::sampling(compiled_bio, data = data_seq, iter = 2000, warmup = 1000, chains = 2, refresh = 0, control = list(adapt_delta = 0.95))
  p_seq   <- apply(plogis(rstan::extract(fit_seq)$alpha), 2, mean)
  b_seq   <- apply(rstan::extract(fit_seq)$beta_cov, 2, mean)
  
  sweep_results <- rbind(sweep_results, data.frame(
    Alpha_Power    = ap,
    Prior_Discount = paste0(round((1 - sqrt(ap)) * 100), "%"),
    Dose_RMSE      = round(rmse(p_seq, true_skeleton), 4),
    Biomarker_AUC  = round(calc_feature_auc(abs(b_seq), causal_mask), 4)
  ))
}

print(sweep_results)