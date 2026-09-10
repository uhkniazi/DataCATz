# ==============================================================================
# File: 06_BIOTITE_PARALLEL_RANKING_SIMULATION.R
# Auth: u.niazi@soton.ac.uk
# Date: 10/09/2026
# Desc: Parallelized Monte Carlo simulation harness evaluating model discrimination,
#       sparsity diagnostic, rank recovery, and Top-K Recall (Top-5, Top-10, Top-20)
#       for high-dimensional Phase I trial data (N = 40, p = 100).
# ==============================================================================

library(rstan)
library(parallel)
library(glmnet)
library(randomForest)
library(xgboost)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. EVALUATION METRICS & DATA GENERATIVE PROCESS (DGP)
# ==============================================================================

#' Calculate Feature Ranking Area Under the ROC Curve (AUC)
#' Converts variable importance scores or coefficient magnitudes into a global AUC
#' evaluating how effectively the model ranks true causal features above noise features.
#'
#' @param scores Numeric vector of feature scores/coefficients (length p = 100).
#' @param true_is_causal Binary indicator vector (1 for true causal signal, 0 for noise).
#' @return Scalar AUC value between 0.0 and 1.0.
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1)
  n0 <- sum(true_is_causal == 0)
  r  <- rank(scores) # Rank features from lowest (1) to highest (p)
  auc <- (sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}

#' Calculate Top-K Recall Fraction
#' Evaluates actionable translational utility: what fraction of the ground-truth
#' causal features (X_1 to X_5) appear within the model's top K ranked candidates.
#'
#' @param beta_scores Numeric vector of feature scores/coefficients.
#' @param true_causal_idx Vector of ground-truth causal feature indices (default 1:5).
#' @param k Integer cutoff for the top feature panel size (e.g., 5, 10, 20).
#' @return Proportion of true causal features captured in the top K list (0.0 to 1.0).
calc_topk_recall <- function(beta_scores, true_causal_idx = 1:5, k = 10) {
  # Sort feature indices by absolute importance score in descending order
  top_k_indices <- order(abs(beta_scores), decreasing = TRUE)[1:k]
  # Count how many ground-truth causal features are captured in top K
  hits <- sum(top_k_indices %in% true_causal_idx)
  return(hits / length(true_causal_idx))
}

#' Stratified Cross-Validation Fold Allocation
#' Ensures balanced response (DLT) representation across cross-validation folds,
#' preventing fold collapse in small sample size (N = 40) settings.
get_stratified_foldid <- function(y, nfolds = 5) {
  foldid <- numeric(length(y))
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  foldid[idx0] <- sample(rep(1:nfolds, length.out = length(idx0)))
  foldid[idx1] <- sample(rep(1:nfolds, length.out = length(idx1)))
  return(foldid)
}

#' Synthesize Phase I Dose-Escalation & High-Dimensional Biomarker Data
#' Generates N = 40 patient profiles with p = 100 features, where 5 features (X_1 to X_5)
#' possess true causal log-odds effects on Dose-Limiting Toxicity (DLT), while 95 are noise.
generate_centered_binary_scenario_data <- function(n = 40, 
                                                   n_doses = 5, 
                                                   n_features = 100, 
                                                   skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                                   dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  # Logit transform dose skeleton probabilities to construct baseline toxicity intercepts
  alpha_true <- qlogis(skeleton)
  
  # Ensure all dose levels are assigned at least one patient
  repeat {
    dose_raw <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)
    if (length(unique(dose_raw)) == n_doses) break
  }
  dose_fac <- factor(dose_raw, levels = 1:n_doses)
  
  # Ground-truth causal effect sizes for X_1 to X_5
  beta_true <- rep(0.0, n_features)
  names(beta_true) <- paste0("X_", 1:n_features)
  beta_true[1:5] <- c(1.0, 0.8, 0.6, 0.5, 0.4) 
  
  X_mat <- matrix(0, nrow = n, ncol = n_features)
  colnames(X_mat) <- paste0("X_", 1:n_features)
  
  # Feature Structure Setup:
  # X_1: Binary, Dose-Confounded
  g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)
  raw_b1  <- rbinom(n, size = 1, prob = g1_prob[dose_raw])
  X_mat[, 1] <- ifelse(raw_b1 == 1, 0.5, -0.5)
  
  # X_2: Continuous, Dose-Confounded
  X_mat[, 2] <- rnorm(n, mean = as.numeric(scale(dose_raw)), sd = 1.0)
  
  # X_3: Continuous, Unconfounded
  X_mat[, 3] <- rnorm(n, mean = 0, sd = 1.0)
  
  # X_4: Binary, Unconfounded
  raw_b4  <- rbinom(n, size = 1, prob = 0.35)
  X_mat[, 4] <- ifelse(raw_b4 == 1, 0.5, -0.5)
  
  # X_5: Continuous, Unconfounded
  X_mat[, 5] <- rnorm(n, mean = 0, sd = 1.0)
  
  # X_6 to X_100: Ground-Truth Noise Variables (Alternating Continuous and Binary)
  for (j in 6:n_features) {
    if (j %% 2 == 0) {
      X_mat[, j] <- rnorm(n, mean = 0, sd = 1.0)
    } else {
      raw_null <- rbinom(n, size = 1, prob = 0.30)
      X_mat[, j] <- ifelse(raw_null == 1, 0.5, -0.5)
    }
  }
  
  # Compute Bernoulli log-odds probability of DLT response
  eta <- alpha_true[dose_raw] + as.vector(X_mat %*% beta_true)
  p   <- plogis(eta)
  dlt <- rbinom(n, size = 1, prob = p)
  
  df_sim   <- data.frame(patient_id = 1:n, dose = dose_fac, dlt = dlt, X_mat)
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
# 2. SINGLE TRIAL SIMULATION WORKER (WITH SELF-HEALING RETRY LOGIC)
# ==============================================================================

#' Execute Single Simulation Worker Iteration
#' Fits LASSO, Elastic Net, Ridge, Random Forest, XGBoost, Standard Stan, and BioTiTE V3.
#' Wrapped in a tryCatch repeat block to silently auto-heal stochastic fitting errors.
run_single_trial_simulation <- function(sim_id, 
                                        compiled_std, 
                                        compiled_biotite, 
                                        n = 40, 
                                        skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50)) {
  
  options(mc.cores = 1)
  
  repeat {
    res <- tryCatch({
      
      # Sample trial data until DLT counts fall in a non-degenerate range (6 <= DLT <= 34)
      repeat {
        sim_data <- generate_centered_binary_scenario_data(n = n, skeleton = skeleton)
        dlt_sum  <- sum(sim_data$df_sim$dlt)
        if (dlt_sum >= 6 && dlt_sum <= (n - 6)) break
      }
      
      df_sim      <- sim_data$df_sim
      X_design    <- sim_data$X_design
      beta_true   <- sim_data$beta_true
      causal_mask <- ifelse(beta_true > 0, 1, 0)
      
      pf     <- c(rep(0, 5), rep(1, 100)) # Unpenalize 5 dose baseline parameters
      foldid <- get_stratified_foldid(df_sim$dlt, nfolds = 5)
      
      # ------------------------------------------------------------------------
      # A. Sparsity Diagnostic Models (LASSO & Elastic Net)
      # ------------------------------------------------------------------------
      fit_lasso <- suppressWarnings(
        cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 1.0, penalty.factor = pf, foldid = foldid)
      )
      fit_enet <- suppressWarnings(
        cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 0.5, penalty.factor = pf, foldid = foldid)
      )
      
      lasso_selected <- sum(as.matrix(coef(fit_lasso, s = "lambda.min"))[paste0("X_", 1:100), 1] != 0)
      enet_selected  <- sum(as.matrix(coef(fit_enet,  s = "lambda.min"))[paste0("X_", 1:100), 1] != 0)
      
      # ------------------------------------------------------------------------
      # B. Ridge Regression
      # ------------------------------------------------------------------------
      fit_ridge  <- suppressWarnings(
        cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 0.0, penalty.factor = pf, foldid = foldid)
      )
      ridge_beta <- as.matrix(coef(fit_ridge, s = "lambda.min"))[paste0("X_", 1:100), 1]
      
      # ------------------------------------------------------------------------
      # C. Random Forest
      # ------------------------------------------------------------------------
      rf_input_df <- df_sim[, c("dose", paste0("X_", 1:100))]
      fit_rf      <- randomForest::randomForest(x = rf_input_df, y = factor(df_sim$dlt), importance = TRUE, ntree = 300)
      rf_imp_mat  <- randomForest::importance(fit_rf, type = 1)
      rf_beta     <- rf_imp_mat[paste0("X_", 1:100), ncol(rf_imp_mat)]
      
      # ------------------------------------------------------------------------
      # D. XGBoost
      # ------------------------------------------------------------------------
      xgb_x  <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:100))])
      dtrain <- xgboost::xgb.DMatrix(data = xgb_x, label = df_sim$dlt)
      fit_xgb <- xgboost::xgb.train(
        params  = list(objective = "binary:logistic", max_depth = 3, learning_rate = 0.1),
        data    = dtrain, nrounds = 50, verbose = 0
      )
      xgb_imp_df <- xgboost::xgb.importance(model = fit_xgb)
      xgb_beta   <- rep(0, 100)
      names(xgb_beta) <- paste0("X_", 1:100)
      matched_xgb <- intersect(names(xgb_beta), xgb_imp_df$Feature)
      if (length(matched_xgb) > 0) {
        xgb_beta[matched_xgb] <- xgb_imp_df$Gain[match(matched_xgb, xgb_imp_df$Feature)]
      }
      
      # ------------------------------------------------------------------------
      # E. Standard Stan GLM
      # ------------------------------------------------------------------------
      stan_data_std <- list(Ntotal = nrow(X_design), Ncol = ncol(X_design), X = X_design, y = as.array(df_sim$dlt))
      fit_std  <- rstan::sampling(compiled_std, data = stan_data_std, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0)
      ext_std  <- rstan::extract(fit_std)
      std_beta <- apply(ext_std$betas[, 6:105], 2, mean)
      names(std_beta) <- paste0("X_", 1:100)
      
      # ------------------------------------------------------------------------
      # F. BioTiTE V3
      # ------------------------------------------------------------------------
      stan_data_bio <- list(
        Ntotal = nrow(X_design), Ncol = ncol(X_design), X = X_design, y = as.array(df_sim$dlt),
        w = rep(1.0, nrow(df_sim)), prior_means = qlogis(skeleton), prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
      )
      fit_bio  <- rstan::sampling(compiled_biotite, data = stan_data_bio, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0)
      ext_bio  <- rstan::extract(fit_bio)
      bio_beta <- apply(ext_bio$beta_cov, 2, mean)
      names(bio_beta) <- paste0("X_", 1:100)
      
      # ------------------------------------------------------------------------
      # G. Feature Ranks & Metrics Calculation
      # ------------------------------------------------------------------------
      ridge_ranks <- rank(abs(ridge_beta))
      rf_ranks    <- rank(rf_beta)
      xgb_ranks   <- rank(xgb_beta)
      std_ranks   <- rank(abs(std_beta))
      bio_ranks   <- rank(abs(bio_beta))
      
      # Global Area Under the Curve (AUC)
      aucs <- c(
        Ridge   = calc_feature_auc(abs(ridge_beta), causal_mask),
        RanFor  = calc_feature_auc(rf_beta, causal_mask),
        XGB     = calc_feature_auc(xgb_beta, causal_mask),
        StdStan = calc_feature_auc(abs(std_beta), causal_mask),
        BioTiTE = calc_feature_auc(abs(bio_beta), causal_mask)
      )
      
      # Top-K Recall Metrics
      top5_rec <- c(
        Ridge   = calc_topk_recall(abs(ridge_beta), 1:5, k = 5),
        RanFor  = calc_topk_recall(rf_beta, 1:5, k = 5),
        XGB     = calc_topk_recall(xgb_beta, 1:5, k = 5),
        StdStan = calc_topk_recall(abs(std_beta), 1:5, k = 5),
        BioTiTE = calc_topk_recall(abs(bio_beta), 1:5, k = 5)
      )
      
      top10_rec <- c(
        Ridge   = calc_topk_recall(abs(ridge_beta), 1:5, k = 10),
        RanFor  = calc_topk_recall(rf_beta, 1:5, k = 10),
        XGB     = calc_topk_recall(xgb_beta, 1:5, k = 10),
        StdStan = calc_topk_recall(abs(std_beta), 1:5, k = 10),
        BioTiTE = calc_topk_recall(abs(bio_beta), 1:5, k = 10)
      )
      
      top20_rec <- c(
        Ridge   = calc_topk_recall(abs(ridge_beta), 1:5, k = 20),
        RanFor  = calc_topk_recall(rf_beta, 1:5, k = 20),
        XGB     = calc_topk_recall(xgb_beta, 1:5, k = 20),
        StdStan = calc_topk_recall(abs(std_beta), 1:5, k = 20),
        BioTiTE = calc_topk_recall(abs(bio_beta), 1:5, k = 20)
      )
      
      list(
        sim_id         = sim_id,
        lasso_selected = lasso_selected,
        enet_selected  = enet_selected,
        ridge_ranks    = ridge_ranks[1:5],
        rf_ranks       = rf_ranks[1:5],
        xgb_ranks      = xgb_ranks[1:5],
        std_ranks      = std_ranks[1:5],
        bio_ranks      = bio_ranks[1:5],
        aucs           = aucs,
        top5_rec       = top5_rec,
        top10_rec      = top10_rec,
        top20_rec      = top20_rec
      )
      
    }, error = function(e) {
      # Return NULL on stochastic failure to trigger clean iteration retry
      NULL
    })
    
    if (!is.null(res)) break
  }
  
  return(res)
}

# ==============================================================================
# 3. METRICS AGGREGATION & SUMMARY GENERATION
# ==============================================================================

#' Aggregate Simulation Metrics Across Monte Carlo Runs
#' Compiles feature selection counts, global AUCs, Top-K recall rates, and causal ranks.
summarize_simulation_results <- function(raw_results) {
  n_sims <- length(raw_results)
  
  # 1. Sparsity Diagnostics
  lasso_sel <- mean(sapply(raw_results, function(x) x$lasso_selected))
  enet_sel  <- mean(sapply(raw_results, function(x) x$enet_selected))
  
  # 2. AUC Summaries
  auc_mat <- t(sapply(raw_results, function(x) x$aucs))
  auc_summary <- data.frame(
    Model_Class  = colnames(auc_mat),
    Mean_AUC     = round(colMeans(auc_mat), 4),
    SE_AUC       = round(apply(auc_mat, 2, sd) / sqrt(n_sims), 4),
    Median_AUC   = round(apply(auc_mat, 2, median), 4)
  )
  
  # 3. Top-K Recall Summaries
  top5_mat  <- t(sapply(raw_results, function(x) x$top5_rec))
  top10_mat <- t(sapply(raw_results, function(x) x$top10_rec))
  top20_mat <- t(sapply(raw_results, function(x) x$top20_rec))
  
  recall_summary <- data.frame(
    Model_Class    = colnames(top5_mat),
    Mean_Top5_Rec  = round(colMeans(top5_mat), 4),
    Mean_Top10_Rec = round(colMeans(top10_mat), 4),
    Mean_Top20_Rec = round(colMeans(top20_mat), 4)
  )
  
  # 4. Mean Causal Ranks (X_1 to X_5)
  ridge_r_mat <- t(sapply(raw_results, function(x) x$ridge_ranks))
  rf_r_mat    <- t(sapply(raw_results, function(x) x$rf_ranks))
  xgb_r_mat   <- t(sapply(raw_results, function(x) x$xgb_ranks))
  std_r_mat   <- t(sapply(raw_results, function(x) x$std_ranks))
  bio_r_mat   <- t(sapply(raw_results, function(x) x$bio_ranks))
  
  rank_summary <- data.frame(
    Feature       = paste0("X_", 1:5),
    Ridge_MeanR   = round(colMeans(ridge_r_mat), 1),
    RanFor_MeanR  = round(colMeans(rf_r_mat), 1),
    XGB_MeanR     = round(colMeans(xgb_r_mat), 1),
    StdStan_MeanR = round(colMeans(std_r_mat), 1),
    BioTiTE_MeanR = round(colMeans(bio_r_mat), 1)
  )
  
  return(list(
    Sparsity_Diagnostic = data.frame(
      Model = c("LASSO (alpha=1.0)", "Elastic Net (alpha=0.5)"),
      Mean_Features_Selected = c(lasso_sel, enet_sel)
    ),
    AUC_Summary         = auc_summary,
    TopK_Recall_Summary = recall_summary,
    Causal_Rank_Summary = rank_summary
  ))
}

# ==============================================================================
# 4. PARALLEL EXECUTION HARNESS
# ==============================================================================

#' Run Parallel Simulation Benchmark
#' Manages sequential Stan compilation, multi-core worker cluster creation,
#' workload distribution, and cluster teardown.
run_parallel_ranking_study <- function(N_sims = 100, N_patients = 40) {
  
  skeleton <- c(0.05, 0.10, 0.20, 0.35, 0.50)
  
  cat("Pre-compiling Stan models sequentially on main node...\n")
  rstan_options(auto_write = TRUE)
  compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
  compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")
  
  num_workers <- min(parallel::detectCores() - 2, N_sims)
  cat(sprintf("Launching parallel cluster across %d workers for %d simulations...\n", num_workers, N_sims))
  cl <- parallel::makeCluster(num_workers)
  
  parallel::clusterEvalQ(cl, {
    library(rstan)
    library(glmnet)
    library(randomForest)
    library(xgboost)
    rstan_options(auto_write = TRUE)
  })
  
  parallel::clusterExport(cl, c(
    "calc_feature_auc", "calc_topk_recall", "get_stratified_foldid", 
    "generate_centered_binary_scenario_data", "run_single_trial_simulation", 
    "N_patients", "skeleton", "compiled_std", "compiled_biotite"
  ), envir = environment())
  
  t_start <- Sys.time()
  
  raw_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
    run_single_trial_simulation(
      sim_id           = s,
      compiled_std     = compiled_std,
      compiled_biotite = compiled_biotite,
      n                = N_patients,
      skeleton         = skeleton
    )
  })
  
  parallel::stopCluster(cl)
  cat(sprintf("Completed %d Monte Carlo simulations in %.2f minutes.\n", N_sims, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  
  summary_out <- summarize_simulation_results(raw_results)
  return(summary_out)
}

# ==============================================================================
# 5. EXECUTION ENTRY POINT
# ==============================================================================
set.seed(42)
sim_results <- run_parallel_ranking_study(N_sims = 500, N_patients = 40)

cat("\n======================================================================\n")
cat("SPARSITY DIAGNOSTIC (Mean Features Selected out of 100)\n")
cat("======================================================================\n")
print(sim_results$Sparsity_Diagnostic)

cat("\n======================================================================\n")
cat("MONTE CARLO AUC SUMMARY (Across Repeated Simulations)\n")
cat("======================================================================\n")
print(sim_results$AUC_Summary)

cat("\n======================================================================\n")
cat("TOP-K CAUSAL FEATURE RECALL SUMMARY (Fraction Captured out of 5)\n")
cat("======================================================================\n")
print(sim_results$TopK_Recall_Summary)

cat("\n======================================================================\n")
cat("MEAN CAUSAL FEATURE RANKS (X_1 to X_5; Higher is better)\n")
cat("======================================================================\n")
print(sim_results$Causal_Rank_Summary)