# ==============================================================================
# File: 05_BIOTITE_PARALLEL_SIMULATION_HARNESS.R
# Auth: u.niazi@soton.ac.uk
# Date: 08/09/2026
# Desc: Parallel simulation harness for BioTiTE Stage 2 evaluation with safe
#       GLM coefficient extraction and empirical dose monotonicity checking.
# ==============================================================================

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. GENERATIVE PROCESS & DATA ASSEMBLY FUNCTION
# ==============================================================================
generate_scenario_data <- function(n = 40, 
                                   n_doses = 5, 
                                   n_genes = 5, 
                                   skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                   beta_true = c(1.0, 0.5, 0.0, 0.0, 0.0),
                                   dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10),
                                   g1_prob = c(0.10, 0.20, 0.40, 0.60, 0.80)) {
  
  alpha_true <- qlogis(skeleton)
  dose_raw   <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
  
  G <- matrix(rnorm(n * n_genes), nrow = n, ncol = n_genes)
  G[, 1] <- rbinom(n, size = 1, prob = g1_prob[dose_raw]) # Dose-confounded Gene1
  colnames(G) <- paste0("Gene", 1:n_genes)
  
  eta <- alpha_true[dose_raw] + as.vector(G %*% beta_true)
  p   <- plogis(eta)
  dlt <- rbinom(n, size = 1, prob = p)
  
  # Explicitly lock factor levels 1:5 to prevent dropping empty dose cohorts
  df_sim <- data.frame(
    patient_id = 1:n, 
    dose       = factor(dose_raw, levels = 1:n_doses), 
    G, 
    p_true     = p, 
    dlt        = dlt
  )
  
  X_design <- model.matrix(~ 0 + dose + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, data = df_sim)
  
  return(list(
    df_sim     = df_sim,
    X_design   = X_design,
    p_true     = p,
    alpha_true = alpha_true
  ))
}

# ==============================================================================
# 2. SINGLE TRIAL SIMULATION WORKER FUNCTION
# ==============================================================================
run_single_trial_simulation <- function(sim_id, 
                                        compiled_std, 
                                        compiled_biotite, 
                                        n = 40, 
                                        skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                        beta_true = c(1.0, 0.5, 0.0, 0.0, 0.0)) {
  
  # Force worker process to strictly use 1 core for MCMC sampling
  options(mc.cores = 1)
  
  # Generate synthetic cohort
  sim_data <- generate_scenario_data(n = n, skeleton = skeleton, beta_true = beta_true)
  df_sim   <- sim_data$df_sim
  X_design <- sim_data$X_design
  
  stan_data_std <- list(
    Ntotal = nrow(X_design),
    Ncol   = ncol(X_design),
    X      = X_design,
    y      = as.array(df_sim$dlt)
  )
  
  stan_data_biotite <- list(
    Ntotal      = nrow(X_design),
    Ncol        = ncol(X_design),
    X           = X_design,
    y           = as.array(df_sim$dlt),
    w           = rep(1.0, nrow(df_sim)),
    prior_means = qlogis(skeleton)
  )
  
  # --- TIER 0: GLM (Safe Name-Based Extraction) ---
  fit_glm <- suppressWarnings(
    glm(dlt ~ 0 + dose + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, 
        data = df_sim, family = binomial(link = "logit"))
  )
  glm_coef_mat <- summary(fit_glm)$coefficients
  
  glm_beta   <- rep(NA_real_, 5)
  glm_dose_p <- rep(NA_real_, 5)
  
  dose_target_names <- paste0("dose", 1:5)
  gene_target_names <- paste0("Gene", 1:5)
  
  for (d in 1:5) {
    if (dose_target_names[d] %in% rownames(glm_coef_mat)) {
      glm_dose_p[d] <- plogis(glm_coef_mat[dose_target_names[d], "Estimate"])
    }
  }
  
  for (g in 1:5) {
    if (gene_target_names[g] %in% rownames(glm_coef_mat)) {
      glm_beta[g] <- glm_coef_mat[gene_target_names[g], "Estimate"]
    }
  }
  
  # --- TIER 1: Standard Stan ---
  fit_std <- rstan::sampling(
    compiled_std, data = stan_data_std, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0
  )
  ext_std  <- rstan::extract(fit_std)
  std_beta <- apply(ext_std$betas[, 6:10], 2, mean)
  std_dose_p <- apply(plogis(ext_std$betas[, 1:5]), 2, mean)
  
  # --- TIER 2: BioTiTE V2 ---
  fit_biotite <- rstan::sampling(
    compiled_biotite, data = stan_data_biotite, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0
  )
  ext_bio  <- rstan::extract(fit_biotite)
  bio_beta <- apply(ext_bio$beta_cov, 2, mean)
  bio_dose_p <- apply(plogis(ext_bio$alpha), 2, mean)
  
  # 89% Interval Coverage calculation
  bio_ci_lower <- apply(ext_bio$beta_cov, 2, quantile, probs = 0.055)
  bio_ci_upper <- apply(ext_bio$beta_cov, 2, quantile, probs = 0.945)
  bio_coverage <- (beta_true >= bio_ci_lower) & (beta_true <= bio_ci_upper)
  
  std_ci_lower <- apply(ext_std$betas[, 6:10], 2, quantile, probs = 0.055)
  std_ci_upper <- apply(ext_std$betas[, 6:10], 2, quantile, probs = 0.945)
  std_coverage <- (beta_true >= std_ci_lower) & (beta_true <= std_ci_upper)
  
  # Explicit empirical monotonicity violation checks across all tiers
  glm_p_valid   <- glm_dose_p[!is.na(glm_dose_p)]
  glm_mono_viol <- if (length(glm_p_valid) > 1) any(diff(glm_p_valid) < 0) else NA
  std_mono_viol <- any(diff(std_dose_p) < 0)
  bio_mono_viol <- any(diff(bio_dose_p) < 0)
  
  return(list(
    sim_id        = sim_id,
    glm_beta      = glm_beta,
    std_beta      = std_beta,
    bio_beta      = bio_beta,
    std_coverage  = std_coverage,
    bio_coverage  = bio_coverage,
    glm_dose_p    = glm_dose_p,
    std_dose_p    = std_dose_p,
    bio_dose_p    = bio_dose_p,
    glm_mono_viol = glm_mono_viol,
    std_mono_viol = std_mono_viol,
    bio_mono_viol = bio_mono_viol,
    tau_median    = median(ext_bio$tau)
  ))
}

# ==============================================================================
# 3. STAGE 2 METRICS AGGREGATION FUNCTION
# ==============================================================================
summarize_stage2_metrics <- function(raw_results, beta_true, true_skeleton) {
  N_sims  <- length(raw_results)
  n_genes <- length(beta_true)
  
  # Extract Coefficient Matrices (Sims x Genes)
  glm_betas <- t(sapply(raw_results, function(x) x$glm_beta))
  std_betas <- t(sapply(raw_results, function(x) x$std_beta))
  bio_betas <- t(sapply(raw_results, function(x) x$bio_beta))
  
  std_cov   <- t(sapply(raw_results, function(x) x$std_coverage))
  bio_cov   <- t(sapply(raw_results, function(x) x$bio_coverage))
  
  # Extract Dose Risk Matrices (Sims x Doses)
  glm_doses <- t(sapply(raw_results, function(x) x$glm_dose_p))
  std_doses <- t(sapply(raw_results, function(x) x$std_dose_p))
  bio_doses <- t(sapply(raw_results, function(x) x$bio_dose_p))
  
  # --- Stage 2 Biomarker Diagnostics ---
  calc_metrics <- function(est_matrix, cov_matrix = NULL) {
    bias <- colMeans(est_matrix, na.rm = TRUE) - beta_true
    rmse <- sqrt(colMeans(sweep(est_matrix, 2, beta_true, "-")^2, na.rm = TRUE))
    null_shrinkage <- mean(abs(est_matrix[, 3:5]), na.rm = TRUE)
    coverage <- if (!is.null(cov_matrix)) colMeans(cov_matrix, na.rm = TRUE) else rep(NA, n_genes)
    return(data.frame(Bias = round(bias, 3), RMSE = round(rmse, 3), Coverage89 = round(coverage, 3), NullIndex = round(null_shrinkage, 3)))
  }
  
  summary_glm <- calc_metrics(glm_betas)
  summary_std <- calc_metrics(std_betas, std_cov)
  summary_bio <- calc_metrics(bio_betas, bio_cov)
  
  biomarker_summary <- data.frame(
    True_Effect = beta_true,
    GLM_Bias = summary_glm$Bias, GLM_RMSE = summary_glm$RMSE,
    Std_Bias = summary_std$Bias, Std_RMSE = summary_std$RMSE, Std_Cov89 = summary_std$Coverage89,
    Bio_Bias = summary_bio$Bias, Bio_RMSE = summary_bio$RMSE, Bio_Cov89 = summary_bio$Coverage89
  )
  rownames(biomarker_summary) <- paste0("Gene", 1:n_genes)
  
  # --- Stage 2 Dose Diagnostics ---
  calc_dose_mse <- function(dose_matrix) {
    mean(sweep(dose_matrix, 2, true_skeleton, "-")^2, na.rm = TRUE)
  }
  
  dose_summary <- data.frame(
    Model = c("GLM (Tier 0)", "Standard Stan (Tier 1)", "BioTiTE V2 (Tier 2)"),
    Dose_MSE = c(calc_dose_mse(glm_doses), calc_dose_mse(std_doses), calc_dose_mse(bio_doses)),
    Monotonicity_Violation_Rate = c(
      mean(sapply(raw_results, function(x) x$glm_mono_viol), na.rm = TRUE),
      mean(sapply(raw_results, function(x) x$std_mono_viol), na.rm = TRUE),
      mean(sapply(raw_results, function(x) x$bio_mono_viol), na.rm = TRUE)
    )
  )
  
  tau_vector <- sapply(raw_results, function(x) x$tau_median)
  
  return(list(
    Biomarker_Recovery = biomarker_summary,
    Dose_Recovery      = dose_summary,
    Null_Shrinkage     = c(GLM = summary_glm$NullIndex[1], StdStan = summary_std$NullIndex[1], BioTiTE = summary_bio$NullIndex[1]),
    Tau_Distribution   = summary(tau_vector)
  ))
}

# ==============================================================================
# 4. PARALLEL CLUSTER EXECUTION HARNESS
# ==============================================================================
run_parallel_simulation_study <- function(N_sims = 100, N_patients = 40) {
  
  cat("Compiling Stan models on main node...\n")
  compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
  compiled_biotite <- rstan::stan_model(file = "BioTiTE_V2.stan")
  
  skeleton  <- c(0.05, 0.10, 0.20, 0.35, 0.50)
  beta_true <- c(1.0, 0.5, 0.0, 0.0, 0.0)
  
  num_workers <- min(parallel::detectCores() - 4, N_sims)
  cat(sprintf("Launching parallel cluster across %d workers...\n", num_workers))
  
  cl <- parallel::makeCluster(num_workers)
  
  parallel::clusterEvalQ(cl, { library(rstan) })
  
  parallel::clusterExport(cl, c(
    "generate_scenario_data", "run_single_trial_simulation", 
    "compiled_std", "compiled_biotite", "N_patients", "skeleton", "beta_true"
  ), envir = environment())
  
  cat(sprintf("Executing %d Monte Carlo simulations...\n", N_sims))
  t_start <- Sys.time()
  
  raw_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
    run_single_trial_simulation(
      sim_id           = s, 
      compiled_std     = compiled_std, 
      compiled_biotite = compiled_biotite, 
      n                = N_patients, 
      skeleton         = skeleton, 
      beta_true        = beta_true
    )
  })
  
  parallel::stopCluster(cl)
  cat(sprintf("Completed in %.2f minutes.\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  
  results <- summarize_stage2_metrics(raw_results, beta_true, skeleton)
  return(results)
}

# ==============================================================================
# 5. EXECUTION ENTRY POINT
# ==============================================================================
set.seed(42)
sim_results <- run_parallel_simulation_study(N_sims = 100, N_patients = 40)
print(sim_results$Dose_Recovery)
print(sim_results$Biomarker_Recovery)