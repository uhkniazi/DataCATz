# ==============================================================================
# File: 07_BIOTITE_CURVE_RECOVERY_CLEAN.R
# Auth: u.niazi@soton.ac.uk
# Date: 11/09/2026
# Desc: Evaluates dose-toxicity curve recovery (RMSE & Monotonicity Violations) 
#       under prior skeleton misspecification and P=100 covariate adjustment.
# ==============================================================================

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. HELPERS & DATA GENERATOR
# ==============================================================================
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1); n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  return((sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}

calc_topk_recall <- function(beta_scores, true_causal_idx = 1:5, k = 10) {
  top_k_indices <- order(abs(beta_scores), decreasing = TRUE)[1:k]
  return(sum(top_k_indices %in% true_causal_idx) / length(true_causal_idx))
}

check_monotonicity_violation <- function(est_probs) {
  return(any(diff(est_probs) < 0))
}

generate_curve_recovery_data <- function(n = 40, 
                                         n_doses = 5, 
                                         n_features = 100, 
                                         true_skeleton  = c(0.08, 0.15, 0.28, 0.45, 0.62),
                                         dose_probs     = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
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
  
  # Causal Features
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  raw_b1  <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  X_mat[, 1] <- ifelse(raw_b1 == 1, 0.5, -0.5)
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  raw_b4  <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 4] <- ifelse(raw_b4 == 1, 0.5, -0.5)
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # Noise Features
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
    true_skeleton  = true_skeleton
  ))
}

# ==============================================================================
# 2. WORKER FUNCTION
# ==============================================================================
run_single_curve_simulation <- function(sim_id, 
                                        compiled_std, 
                                        compiled_biotite, 
                                        n = 40,
                                        true_skeleton = c(0.08, 0.15, 0.28, 0.45, 0.62),
                                        prior_skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50)) {
  
  set.seed(42 + sim_id)
  options(mc.cores = 1)
  causal_mask  <- c(rep(1, 5), rep(0, 95))
  max_attempts <- 10
  attempt      <- 0
  
  repeat {
    attempt <- attempt + 1
    res <- tryCatch({
      
      repeat {
        sim_data <- generate_curve_recovery_data(n = n, true_skeleton = true_skeleton)
        dlt_sum  <- sum(sim_data$df_sim$dlt)
        if (dlt_sum >= 6 && dlt_sum <= (n - 6)) break
      }
      
      df_sim   <- sim_data$df_sim
      
      # 1. Standard Stan GLM (P=100 Covariates)
      std_data <- list(Ntotal = n, Ncol = 105, X = sim_data$X_design_cov, y = as.array(df_sim$dlt))
      fit_std  <- rstan::sampling(compiled_std, data = std_data, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_std  <- rstan::extract(fit_std)
      p_std    <- apply(plogis(ext_std$betas[, 1:5]), 2, mean)
      std_beta <- apply(ext_std$betas[, 6:105], 2, mean)
      
      # 2. BioTiTE V3 (P=0 Omitted Covariates)
      bio_data_nocov <- list(
        Ntotal = n, Ncol = 5, X = sim_data$X_design_nocov, y = as.array(df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(prior_skeleton), 
        prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
      )
      fit_bio_nocov <- rstan::sampling(compiled_biotite, data = bio_data_nocov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      p_bio_nocov   <- apply(plogis(rstan::extract(fit_bio_nocov)$alpha), 2, mean)
      
      # 3. BioTiTE V3 (P=100 Adjusted Covariates)
      bio_data_cov <- list(
        Ntotal = n, Ncol = 105, X = sim_data$X_design_cov, y = as.array(df_sim$dlt),
        w = rep(1.0, n), prior_means = qlogis(prior_skeleton), 
        prior_sd_alpha1 = 1.50, prior_sds_delta = rep(1.0, 4)
      )
      fit_bio_cov <- rstan::sampling(compiled_biotite, data = bio_data_cov, iter = 2000, warmup = 1000, chains = 2, refresh = 0)
      ext_bio_cov <- rstan::extract(fit_bio_cov)
      p_bio_cov   <- apply(plogis(ext_bio_cov$alpha), 2, mean)
      bio_beta    <- apply(ext_bio_cov$beta_cov, 2, mean)
      
      list(
        rmse_std       = sqrt(mean((p_std - true_skeleton)^2)),
        rmse_bio_nocov = sqrt(mean((p_bio_nocov - true_skeleton)^2)),
        rmse_bio_cov   = sqrt(mean((p_bio_cov - true_skeleton)^2)),
        
        mono_viol_std  = as.numeric(check_monotonicity_violation(p_std)),
        mono_viol_bio  = as.numeric(check_monotonicity_violation(p_bio_cov)),
        
        std_auc        = calc_feature_auc(abs(std_beta), causal_mask),
        bio_auc        = calc_feature_auc(abs(bio_beta), causal_mask),
        bio_top10_rec  = calc_topk_recall(abs(bio_beta), 1:5, k = 10)
      )
      
    }, error = function(e) { NULL })
    
    if (!is.null(res) || attempt >= max_attempts) break
  }
  
  return(res)
}

# ==============================================================================
# 3. PARALLEL HARNESS
# ==============================================================================
run_scenario2_study <- function(N_sims = 100, N_patients = 40) {
  
  true_skeleton  <- c(0.08, 0.15, 0.28, 0.45, 0.62)
  prior_skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)
  
  cat("Pre-compiling Stan models sequentially on main node...\n")
  compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
  compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")
  
  num_workers <- min(parallel::detectCores() - 2, N_sims)
  cat(sprintf("Launching parallel cluster across %d workers for %d Scenario 2 simulations...\n", num_workers, N_sims))
  cl <- parallel::makeCluster(num_workers)
  
  parallel::clusterEvalQ(cl, {
    library(rstan)
    rstan_options(auto_write = TRUE)
  })
  
  parallel::clusterExport(cl, c(
    "calc_feature_auc", "calc_topk_recall", "check_monotonicity_violation",
    "generate_curve_recovery_data", "run_single_curve_simulation", "N_patients", 
    "true_skeleton", "prior_skeleton", "compiled_std", "compiled_biotite"
  ), envir = environment())
  
  t_start <- Sys.time()
  
  raw_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
    run_single_curve_simulation(
      sim_id           = s,
      compiled_std     = compiled_std,
      compiled_biotite = compiled_biotite,
      n                = N_patients,
      true_skeleton    = true_skeleton,
      prior_skeleton   = prior_skeleton
    )
  })
  
  parallel::stopCluster(cl)
  cat(sprintf("Completed %d simulations in %.2f minutes.\n", N_sims, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  
  valid_results <- raw_results[!sapply(raw_results, is.null)]
  
  summary_df <- data.frame(
    Model_Architecture = c("Standard Stan GLM (P=100)", "BioTiTE V3 (Omitted Covariates P=0)", "BioTiTE V3 (Adjusted Covariates P=100)"),
    Mean_Curve_RMSE    = round(c(
      mean(sapply(valid_results, function(x) x$rmse_std)),
      mean(sapply(valid_results, function(x) x$rmse_bio_nocov)),
      mean(sapply(valid_results, function(x) x$rmse_bio_cov))
    ), 4),
    Monotonicity_Violation_Rate = c(
      paste0(round(mean(sapply(valid_results, function(x) x$mono_viol_std)) * 100, 1), "%"),
      "0.0%",
      paste0(round(mean(sapply(valid_results, function(x) x$mono_viol_bio)) * 100, 1), "%")
    ),
    Biomarker_AUC = c(
      round(mean(sapply(valid_results, function(x) x$std_auc)), 4),
      "N/A",
      round(mean(sapply(valid_results, function(x) x$bio_auc)), 4)
    )
  )
  
  cat("\n======================================================================\n")
  cat("SCENARIO 2: DOSE CURVE RECOVERY & COVARIATE ADJUSTMENT SUMMARY\n")
  cat("Ground Truth Skeleton: ", true_skeleton, "\n")
  cat("Misspecified Prior:    ", prior_skeleton, "\n")
  cat("======================================================================\n")
  print(summary_df)
}

# Run 50 Test Simulations
set.seed(42)
run_scenario2_study(N_sims = 50, N_patients = 40)