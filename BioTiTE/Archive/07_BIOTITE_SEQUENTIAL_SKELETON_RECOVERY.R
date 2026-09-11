# ==============================================================================
# File: 07_BIOTITE_SEQUENTIAL_SKELETON_RECOVERY.R
# Auth: u.niazi@soton.ac.uk
# Date: 10/09/2026
# Desc: Evaluates sequential trial-to-trial dose skeleton learning (Trial 1 -> Trial 2)
#       under misspecified prior skeletons with and without p=100 covariate adjustment.
# ==============================================================================

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. HELPERS & DATA GENERATOR
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
  
  # Feature generation
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  raw_b1  <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  X_mat[, 1] <- ifelse(raw_b1 == 1, 0.5, -0.5)
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  raw_b4  <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 4] <- ifelse(raw_b4 == 1, 0.5, -0.5)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
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
                                      true_skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                      misspecified_skeleton = c(0.02, 0.05, 0.10, 0.20, 0.35)) {
  
  options(mc.cores = 1)
  causal_mask <- c(rep(1, 5), rep(0, 95))
  
  repeat {
    res <- tryCatch({
      
      # ------------------------------------------------------------------------
      # TRIAL 1 EXECUTION (Misspecified Prior Input)
      # ------------------------------------------------------------------------
      repeat {
        t1_data <- generate_scenario2_trial_data(n = n, true_skeleton = true_skeleton)
        if (sum(t1_data$df_sim$dlt) >= 6 && sum(t1_data$df_sim$dlt) <= (n - 6)) break
      }
      
      # Fit BioTiTE T1 (With Covariates)
      bio_t1_data_cov <- list(
        Ntotal = n, Ncol = 105, X = t1_data$X_design_cov, y = as.array(t1_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(misspecified_skeleton), 
        prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
      )
      fit_bio_t1_cov <- rstan::sampling(compiled_biotite, data = bio_t1_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_bio_t1_cov <- rstan::extract(fit_bio_t1_cov)
      
      # Extract Trial 1 Posterior Skeleton (Logit Scale)
      t1_bio_alpha_means <- apply(ext_bio_t1_cov$alpha, 2, mean)
      t1_bio_alpha_sds   <- apply(ext_bio_t1_cov$alpha, 2, sd)
      p_bio_t1_cov       <- plogis(t1_bio_alpha_means)
      
      # Fit BioTiTE T1 (No Covariates)
      bio_t1_data_nocov <- list(
        Ntotal = n, Ncol = 5, X = t1_data$X_design_nocov, y = as.array(t1_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(misspecified_skeleton), 
        prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
      )
      fit_bio_t1_nocov <- rstan::sampling(compiled_biotite, data = bio_t1_data_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_bio_t1_nocov   <- apply(plogis(rstan::extract(fit_bio_t1_nocov)$alpha), 2, mean)
      
      # Fit Standard Stan T1 (With Covariates)
      std_t1_data_cov <- list(Ntotal = n, Ncol = 105, X = t1_data$X_design_cov, y = as.array(t1_data$df_sim$dlt))
      fit_std_t1_cov  <- rstan::sampling(compiled_std, data = std_t1_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_std_t1_cov    <- apply(plogis(rstan::extract(fit_std_t1_cov)$betas[, 1:5]), 2, mean)
      
      # ------------------------------------------------------------------------
      # TRIAL 2 EXECUTION (Inject Trial 1 Posterior as Trial 2 Prior)
      # ------------------------------------------------------------------------
      repeat {
        t2_data <- generate_scenario2_trial_data(n = n, true_skeleton = true_skeleton)
        if (sum(t2_data$df_sim$dlt) >= 6 && sum(t2_data$df_sim$dlt) <= (n - 6)) break
      }
      
      # BioTiTE T2: Updated Prior Means/SDs from T1
      t2_prior_means_delta <- diff(t1_bio_alpha_means)
      t2_prior_sds_delta   <- sqrt(t1_bio_alpha_sds[1:4]^2 + t1_bio_alpha_sds[2:5]^2)
      
      bio_t2_data_cov <- list(
        Ntotal = n, Ncol = 105, X = t2_data$X_design_cov, y = as.array(t2_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = t1_bio_alpha_means, 
        prior_sd_alpha1 = t1_bio_alpha_sds[1], prior_sds_delta = t2_prior_sds_delta
      )
      fit_bio_t2_cov <- rstan::sampling(compiled_biotite, data = bio_t2_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_bio_t2_cov <- rstan::extract(fit_bio_t2_cov)
      
      p_bio_t2_cov   <- apply(plogis(ext_bio_t2_cov$alpha), 2, mean)
      bio_t2_beta    <- apply(ext_bio_t2_cov$beta_cov, 2, mean)
      
      # Fit BioTiTE T2 (No Covariates)
      bio_t2_data_nocov <- list(
        Ntotal = n, Ncol = 5, X = t2_data$X_design_nocov, y = as.array(t2_data$df_sim$dlt),
        w = rep(1.0, n), prior_means = t1_bio_alpha_means, 
        prior_sd_alpha1 = t1_bio_alpha_sds[1], prior_sds_delta = t2_prior_sds_delta
      )
      fit_bio_t2_nocov <- rstan::sampling(compiled_biotite, data = bio_t2_data_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_bio_t2_nocov   <- apply(plogis(rstan::extract(fit_bio_t2_nocov)$alpha), 2, mean)
      
      # Standard Stan T2 (With Covariates)
      std_t2_data_cov <- list(Ntotal = n, Ncol = 105, X = t2_data$X_design_cov, y = as.array(t2_data$df_sim$dlt))
      fit_std_t2_cov  <- rstan::sampling(compiled_std, data = std_t2_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_std_t2_cov    <- apply(plogis(rstan::extract(fit_std_t2_cov)$betas[, 1:5]), 2, mean)
      
      # ------------------------------------------------------------------------
      # METRICS COMPILATION
      # ------------------------------------------------------------------------
      rmse_bio_t1_cov   <- sqrt(mean((p_bio_t1_cov - true_skeleton)^2))
      rmse_bio_t1_nocov <- sqrt(mean((p_bio_t1_nocov - true_skeleton)^2))
      rmse_std_t1_cov   <- sqrt(mean((p_std_t1_cov - true_skeleton)^2))
      
      rmse_bio_t2_cov   <- sqrt(mean((p_bio_t2_cov - true_skeleton)^2))
      rmse_bio_t2_nocov <- sqrt(mean((p_bio_t2_nocov - true_skeleton)^2))
      rmse_std_t2_cov   <- sqrt(mean((p_std_t2_cov - true_skeleton)^2))
      
      t2_auc      <- calc_feature_auc(abs(bio_t2_beta), causal_mask)
      t2_top10_rec <- calc_topk_recall(abs(bio_t2_beta), 1:5, k = 10)
      
      list(
        rmse_bio_t1_cov   = rmse_bio_t1_cov,
        rmse_bio_t1_nocov = rmse_bio_t1_nocov,
        rmse_std_t1_cov   = rmse_std_t1_cov,
        rmse_bio_t2_cov   = rmse_bio_t2_cov,
        rmse_bio_t2_nocov = rmse_bio_t2_nocov,
        rmse_std_t2_cov   = rmse_std_t2_cov,
        t2_auc            = t2_auc,
        t2_top10_rec      = t2_top10_rec
      )
      
    }, error = function(e) { NULL })
    
    if (!is.null(res)) break
  }
  
  return(res)
}

# ==============================================================================
# 3. PARALLEL HARNESS & EXECUTION ENTRY POINT
# ==============================================================================
run_scenario2_study <- function(N_sims = 100, N_patients = 40) {
  
  true_skeleton         <- c(0.05, 0.10, 0.20, 0.35, 0.50)
  misspecified_skeleton <- c(0.02, 0.05, 0.10, 0.20, 0.35)
  
  cat("Pre-compiling Stan models sequentially on main node...\n")
  compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
  compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")
  
  num_workers <- min(parallel::detectCores() - 2, N_sims)
  cat(sprintf("Launching cluster across %d workers for %d Scenario 2 simulations...\n", num_workers, N_sims))
  cl <- parallel::makeCluster(num_workers)
  
  parallel::clusterEvalQ(cl, {
    library(rstan)
    rstan_options(auto_write = TRUE)
  })
  
  parallel::clusterExport(cl, c(
    "calc_feature_auc", "calc_topk_recall", "generate_scenario2_trial_data", 
    "run_sequential_simulation", "N_patients", "true_skeleton", 
    "misspecified_skeleton", "compiled_std", "compiled_biotite"
  ), envir = environment())
  
  t_start <- Sys.time()
  
  raw_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
    run_sequential_simulation(
      sim_id                = s,
      compiled_std          = compiled_std,
      compiled_biotite      = compiled_biotite,
      n                     = N_patients,
      true_skeleton         = true_skeleton,
      misspecified_skeleton = misspecified_skeleton
    )
  })
  
  parallel::stopCluster(cl)
  cat(sprintf("Completed %d simulations in %.2f minutes.\n", N_sims, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  
  # Aggregate RMSE Results
  res_df <- data.frame(
    Architecture = c("Std Stan (Covariates P=100)", "BioTiTE V3 (Omitted Covariates P=0)", "BioTiTE V3 (Adjusted Covariates P=100)"),
    Trial1_RMSE  = c(
      mean(sapply(raw_results, function(x) x$rmse_std_t1_cov)),
      mean(sapply(raw_results, function(x) x$rmse_bio_t1_nocov)),
      mean(sapply(raw_results, function(x) x$rmse_bio_t1_cov))
    ),
    Trial2_RMSE  = c(
      mean(sapply(raw_results, function(x) x$rmse_std_t2_cov)),
      mean(sapply(raw_results, function(x) x$rmse_bio_t2_nocov)),
      mean(sapply(raw_results, function(x) x$rmse_bio_t2_cov))
    )
  )
  res_df$RMSE_Reduction_Pct <- round(((res_df$Trial1_RMSE - res_df$Trial2_RMSE) / res_df$Trial1_RMSE) * 100, 2)
  res_df$Trial1_RMSE <- round(res_df$Trial1_RMSE, 4)
  res_df$Trial2_RMSE <- round(res_df$Trial2_RMSE, 4)
  
  t2_auc_mean   <- round(mean(sapply(raw_results, function(x) x$t2_auc)), 4)
  t2_rec10_mean <- round(mean(sapply(raw_results, function(x) x$t2_top10_rec)), 4)
  
  cat("\n======================================================================\n")
  cat("SCENARIO 2: DOSE SKELETON RECOVERY & SEQUENTIAL LEARNING SUMMARY\n")
  cat("True Skeleton:          ", true_skeleton, "\n")
  cat("Initial Input Skeleton: ", misspecified_skeleton, "\n")
  cat("======================================================================\n")
  print(res_df)
  
  cat(sprintf("\nTrial 2 Biomarker Screening Check (BioTiTE Covariates P=100):\n"))
  cat(sprintf("  Mean AUC:          %.4f\n", t2_auc_mean))
  cat(sprintf("  Mean Top-10 Recall: %.4f\n", t2_rec10_mean))
}

# Run 50 Test Simulations
set.seed(42)
run_scenario2_study(N_sims = 50, N_patients = 40)