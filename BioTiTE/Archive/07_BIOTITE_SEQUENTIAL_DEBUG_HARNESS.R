# ==============================================================================
# File: 07_BIOTITE_SEQUENTIAL_DEBUG_HARNESS.R
# Auth: u.niazi@soton.ac.uk
# Date: 10/09/2026
# Desc: Debug harness for Scenario 2 sequential trial transfer (Trial 1 -> Trial 2)
#       with power-prior variance discount (alpha_power = 0.5), joint covariance 
#       delta extraction, Jensen-corrected posterior expectations, and sub-scenarios 2A, 2B, 2C.
# ==============================================================================

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. EVALUATION METRICS & DATA GENERATION
# ==============================================================================
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1)
  n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  return((sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}

calc_topk_recall <- function(beta_scores, true_causal_idx = 1:5, k = 10) {
  top_k_indices <- order(abs(beta_scores), decreasing = TRUE)[1:k]
  return(sum(top_k_indices %in% true_causal_idx) / length(true_causal_idx))
}

generate_scenario2_trial_data <- function(n = 40, 
                                          n_doses = 5, 
                                          n_features = 100, 
                                          true_skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                          dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  alpha_true <- qlogis(true_skeleton)
  
  repeat {
    dose_raw <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
    if (length(unique(dose_raw)) == n_doses) break
  }
  dose_fac <- factor(dose_raw, levels = 1:n_doses)
  
  beta_true <- rep(0.0, n_features)
  names(beta_true) <- paste0("X_", 1:n_features)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) 
  
  X_mat <- matrix(0, nrow = n, ncol = n_features)
  colnames(X_mat) <- paste0("X_", 1:n_features)
  
  # Feature 1: Binary & Dose-Confounded
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
  
  # Features 6-100: Noise
  for (j in 6:n_features) {
    if (j %% 2 == 0) {
      X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
    } else {
      raw_null <- rbinom(n, size = 1, prob = 0.30)
      X_mat[, j] <- ifelse(raw_null == 1, 0.5, -0.5)
    }
  }
  
  eta <- alpha_true[dose_raw] + as.vector(X_mat %*% beta_true)
  p   <- plogis(eta)
  dlt <- rbinom(n, size = 1, prob = p)
  
  df_sim <- data.frame(patient_id = 1:n, dose = dose_fac, dlt = dlt, X_mat)
  X_design_cov   <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:n_features))])
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
# 2. SEQUENTIAL WORKER (TRIAL 1 -> TRIAL 2)
# ==============================================================================
run_sequential_simulation <- function(sim_id, 
                                      compiled_std, 
                                      compiled_biotite, 
                                      n = 40,
                                      true_skeleton_t1 = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                      true_skeleton_t2 = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                      misspecified_skeleton = c(0.02, 0.05, 0.10, 0.20, 0.35),
                                      alpha_power = 0.5) {
  
  # Per-worker seed for reproducible parallel streams
  set.seed(42 + sim_id)
  options(mc.cores = 1)
  causal_mask  <- c(rep(1, 5), rep(0, 95))
  max_attempts <- 10
  attempt      <- 0
  
  repeat {
    attempt <- attempt + 1
    res <- tryCatch({
      
      # ------------------------------------------------------------------------
      # TRIAL 1: Misspecified Prior Input
      # ------------------------------------------------------------------------
      repeat {
        t1_data <- generate_scenario2_trial_data(n = n, true_skeleton = true_skeleton_t1)
        if (sum(t1_data$df_sim$dlt) >= 6 && sum(t1_data$df_sim$dlt) <= (n - 6)) break
      }
      
      # BioTiTE T1 (With Covariates, Loosened Prior SDs)
      bio_t1_data_cov <- list(
        Ntotal = n, Ncol = 105, X = t1_data$X_design_cov, y = as.array(t1_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(misspecified_skeleton), 
        prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
      )
      fit_bio_t1_cov <- rstan::sampling(compiled_biotite, data = bio_t1_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_bio_t1_cov <- rstan::extract(fit_bio_t1_cov)
      
      # FIX 1: Correct Jensen Expectation for Probability Curve
      p_bio_t1_cov   <- apply(plogis(ext_bio_t1_cov$alpha), 2, mean)
      t1_bio_beta    <- apply(ext_bio_t1_cov$beta_cov, 2, mean)
      t1_auc         <- calc_feature_auc(abs(t1_bio_beta), causal_mask)
      t1_top10_rec   <- calc_topk_recall(abs(t1_bio_beta), 1:5, k = 10)
      
      # BioTiTE T1 (No Covariates)
      bio_t1_data_nocov <- list(
        Ntotal = n, Ncol = 5, X = t1_data$X_design_nocov, y = as.array(t1_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(misspecified_skeleton), 
        prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
      )
      fit_bio_t1_nocov <- rstan::sampling(compiled_biotite, data = bio_t1_data_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_bio_t1_nocov   <- apply(plogis(rstan::extract(fit_bio_t1_nocov)$alpha), 2, mean)
      
      # Standard Stan T1
      std_t1_data_cov <- list(Ntotal = n, Ncol = 105, X = t1_data$X_design_cov, y = as.array(t1_data$df_sim$dlt))
      fit_std_t1_cov  <- rstan::sampling(compiled_std, data = std_t1_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_std_t1_cov    <- apply(plogis(rstan::extract(fit_std_t1_cov)$betas[, 1:5]), 2, mean)
      
      # ------------------------------------------------------------------------
      # FIX 2 & 3: COVARIANCE-AWARE DELTA EXTRACTION & POWER-PRIOR SCALING
      # ------------------------------------------------------------------------
      t1_alpha_samples <- ext_bio_t1_cov$alpha                              # S x 5
      t1_delta_samples <- t1_alpha_samples[, 2:5] - t1_alpha_samples[, 1:4] # Direct joint delta samples
      
      t1_alpha1_mean <- mean(t1_alpha_samples[, 1])
      t1_alpha1_sd   <- sd(t1_alpha_samples[, 1])
      
      t1_delta_means <- colMeans(t1_delta_samples)
      t1_delta_sds   <- apply(t1_delta_samples, 2, sd)
      
      # Apply Power-Prior Variance Inflation (alpha_power = 0.5)
      t2_prior_sd_alpha1   <- t1_alpha1_sd / sqrt(alpha_power)
      t2_prior_sds_delta   <- t1_delta_sds / sqrt(alpha_power)
      t2_prior_means_alpha <- c(t1_alpha1_mean, t1_alpha1_mean + cumsum(t1_delta_means))
      
      # ------------------------------------------------------------------------
      # TRIAL 2: Sequential Transfer Execution
      # ------------------------------------------------------------------------
      repeat {
        t2_data <- generate_scenario2_trial_data(n = n, true_skeleton = true_skeleton_t2)
        if (sum(t2_data$df_sim$dlt) >= 6 && sum(t2_data$df_sim$dlt) <= (n - 6)) break
      }
      
      # BioTiTE T2 (With Covariates)
      bio_t2_data_cov <- list(
        Ntotal = n, Ncol = 105, X = t2_data$X_design_cov, y = as.array(t2_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = t2_prior_means_alpha, 
        prior_sd_alpha1 = t2_prior_sd_alpha1, prior_sds_delta = t2_prior_sds_delta
      )
      fit_bio_t2_cov <- rstan::sampling(compiled_biotite, data = bio_t2_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_bio_t2_cov <- rstan::extract(fit_bio_t2_cov)
      
      p_bio_t2_cov   <- apply(plogis(ext_bio_t2_cov$alpha), 2, mean)
      bio_t2_beta    <- apply(ext_bio_t2_cov$beta_cov, 2, mean)
      t2_auc         <- calc_feature_auc(abs(bio_t2_beta), causal_mask)
      t2_top10_rec   <- calc_topk_recall(abs(bio_t2_beta), 1:5, k = 10)
      
      # BioTiTE T2 (No Covariates)
      bio_t2_data_nocov <- list(
        Ntotal = n, Ncol = 5, X = t2_data$X_design_nocov, y = as.array(t2_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = t2_prior_means_alpha, 
        prior_sd_alpha1 = t2_prior_sd_alpha1, prior_sds_delta = t2_prior_sds_delta
      )
      fit_bio_t2_nocov <- rstan::sampling(compiled_biotite, data = bio_t2_data_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_bio_t2_nocov   <- apply(plogis(rstan::extract(fit_bio_t2_nocov)$alpha), 2, mean)
      
      # Standard Stan T2
      std_t2_data_cov <- list(Ntotal = n, Ncol = 105, X = t2_data$X_design_cov, y = as.array(t2_data$df_sim$dlt))
      fit_std_t2_cov  <- rstan::sampling(compiled_std, data = std_t2_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_std_t2_cov    <- apply(plogis(rstan::extract(fit_std_t2_cov)$betas[, 1:5]), 2, mean)
      
      # ------------------------------------------------------------------------
      # METRICS COMPILATION
      # ------------------------------------------------------------------------
      list(
        rmse_bio_t1_cov   = sqrt(mean((p_bio_t1_cov - true_skeleton_t1)^2)),
        rmse_bio_t1_nocov = sqrt(mean((p_bio_t1_nocov - true_skeleton_t1)^2)),
        rmse_std_t1_cov   = sqrt(mean((p_std_t1_cov - true_skeleton_t1)^2)),
        
        rmse_bio_t2_cov   = sqrt(mean((p_bio_t2_cov - true_skeleton_t2)^2)),
        rmse_bio_t2_nocov = sqrt(mean((p_bio_t2_nocov - true_skeleton_t2)^2)),
        rmse_std_t2_cov   = sqrt(mean((p_std_t2_cov - true_skeleton_t2)^2)),
        
        t1_auc            = t1_auc,
        t1_top10_rec      = t1_top10_rec,
        t2_auc            = t2_auc,
        t2_top10_rec      = t2_top10_rec
      )
      
    }, error = function(e) { NULL })
    
    # FIX 7: Graceful Retry Cap
    if (!is.null(res) || attempt >= max_attempts) break
  }
  
  return(res)
}

# ==============================================================================
# 3. PARALLEL HARNESS & SCENARIO CONTROLLER
# ==============================================================================
run_scenario2_subscenario <- function(sub_code, true_t1, true_t2, misspec, N_sims = 10, compiled_std, compiled_biotite) {
  
  cat(sprintf("\n--- Running Sub-Scenario %s (%d Sims) ---\n", sub_code, N_sims))
  
  num_workers <- min(parallel::detectCores() - 2, N_sims)
  cl <- parallel::makeCluster(num_workers)
  
  parallel::clusterEvalQ(cl, {
    library(rstan)
    rstan_options(auto_write = TRUE)
  })
  
  parallel::clusterExport(cl, c(
    "calc_feature_auc", "calc_topk_recall", "generate_scenario2_trial_data", 
    "run_sequential_simulation", "true_t1", "true_t2", "misspec", 
    "compiled_std", "compiled_biotite"
  ), envir = environment())
  
  raw_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
    run_sequential_simulation(
      sim_id                = s,
      compiled_std          = compiled_std,
      compiled_biotite      = compiled_biotite,
      n                     = 40,
      true_skeleton_t1      = true_t1,
      true_skeleton_t2      = true_t2,
      misspecified_skeleton = misspec,
      alpha_power           = 0.5
    )
  })
  
  parallel::stopCluster(cl)
  
  # Filter out NULL workers
  valid_results <- raw_results[!sapply(raw_results, is.null)]
  
  res_df <- data.frame(
    SubScenario  = sub_code,
    Architecture = c("Std Stan (Covariates P=100)", "BioTiTE V3 (Omitted Covariates P=0)", "BioTiTE V3 (Adjusted Covariates P=100)"),
    Trial1_RMSE  = round(c(
      mean(sapply(valid_results, function(x) x$rmse_std_t1_cov)),
      mean(sapply(valid_results, function(x) x$rmse_bio_t1_nocov)),
      mean(sapply(valid_results, function(x) x$rmse_bio_t1_cov))
    ), 4),
    Trial2_RMSE  = round(c(
      mean(sapply(valid_results, function(x) x$rmse_std_t2_cov)),
      mean(sapply(valid_results, function(x) x$rmse_bio_t2_nocov)),
      mean(sapply(valid_results, function(x) x$rmse_bio_t2_cov))
    ), 4)
  )
  res_df$RMSE_Reduction_Pct <- round(((res_df$Trial1_RMSE - res_df$Trial2_RMSE) / res_df$Trial1_RMSE) * 100, 2)
  
  bio_cov_res <- valid_results[[1]]
  t1_auc      <- round(mean(sapply(valid_results, function(x) x$t1_auc)), 4)
  t2_auc      <- round(mean(sapply(valid_results, function(x) x$t2_auc)), 4)
  t1_rec10    <- round(mean(sapply(valid_results, function(x) x$t1_top10_rec)), 4)
  t2_rec10    <- round(mean(sapply(valid_results, function(x) x$t2_top10_rec)), 4)
  
  return(list(
    summary_table = res_df,
    biomarker_check = data.frame(
      SubScenario = sub_code,
      T1_AUC = t1_auc, T2_AUC = t2_auc,
      T1_Top10_Rec = t1_rec10, T2_Top10_Rec = t2_rec10
    )
  ))
}

# ==============================================================================
# 4. EXECUTION ENTRY POINT (10 SIM DEBUG RUN)
# ==============================================================================
cat("Pre-compiling Stan models sequentially on main node...\n")
compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")

misspecified_skeleton <- c(0.02, 0.05, 0.10, 0.20, 0.35)

# Sub-Scenario 2A: Identical Populations (T1 == T2)
true_t1_2a <- c(0.05, 0.10, 0.20, 0.35, 0.50)
true_t2_2a <- c(0.05, 0.10, 0.20, 0.35, 0.50)

# Sub-Scenario 2B: Moderately Shifted Curve (T2 slightly higher toxicity)
true_t1_2b <- c(0.05, 0.10, 0.20, 0.35, 0.50)
true_t2_2b <- c(0.08, 0.15, 0.25, 0.40, 0.55)

# Sub-Scenario 2C: Substantially Different Population (T2 toxic drift)
true_t1_2c <- c(0.05, 0.10, 0.20, 0.35, 0.50)
true_t2_2c <- c(0.10, 0.20, 0.35, 0.50, 0.70)

res_2a <- run_scenario2_subscenario("2A (Identical)", true_t1_2a, true_t2_2a, misspecified_skeleton, N_sims = 10, compiled_std, compiled_biotite)
res_2b <- run_scenario2_subscenario("2B (Shifted)",   true_t1_2b, true_t2_2b, misspecified_skeleton, N_sims = 10, compiled_std, compiled_biotite)
res_2c <- run_scenario2_subscenario("2C (Drifted)",   true_t1_2c, true_t2_2c, misspecified_skeleton, N_sims = 10, compiled_std, compiled_biotite)

cat("\n======================================================================\n")
cat("SCENARIO 2 DEBUG SUMMARY: DOSE CURVE RECOVERY (RMSE)\n")
cat("======================================================================\n")
print(rbind(res_2a$summary_table, res_2b$summary_table, res_2c$summary_table))

cat("\n======================================================================\n")
cat("SCENARIO 2 DEBUG SUMMARY: BIOMARKER SCREENING RETENTION CHECK\n")
cat("======================================================================\n")
print(rbind(res_2a$biomarker_check, res_2b$biomarker_check, res_2c$biomarker_check))