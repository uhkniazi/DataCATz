# ==============================================================================
# File: 06_DEBUG_SEQUENTIAL_5SIMS.R
# Desc: Sequential 5-simulation diagnostic run with stratified cv.glmnet foldid
# ==============================================================================

library(rstan)
library(glmnet)
library(randomForest)
library(xgboost)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)
set.seed(42)

# ==============================================================================
# 1. HELPERS & DATA GENERATOR
# ==============================================================================
calc_feature_auc <- function(scores, true_is_causal) {
  n1 <- sum(true_is_causal == 1)
  n0 <- sum(true_is_causal == 0)
  r  <- rank(scores)
  auc <- (sum(r[true_is_causal == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}

generate_centered_binary_scenario_data <- function(n = 40, 
                                                   n_doses = 5, 
                                                   n_features = 100, 
                                                   skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50),
                                                   dose_probs = c(0.30, 0.30, 0.20, 0.10, 0.10)) {
  
  alpha_true <- qlogis(skeleton)
  
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
  X_design <- model.matrix(~ 0 + dose + ., data = df_sim[, c("dose", paste0("X_", 1:n_features))])
  
  return(list(
    df_sim     = df_sim,
    X_design   = X_design,
    beta_true  = beta_true,
    alpha_true = alpha_true,
    skeleton   = skeleton
  ))
}

# Helper to build stratified fold IDs for binary outcomes
get_stratified_foldid <- function(y, nfolds = 5) {
  foldid <- numeric(length(y))
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  foldid[idx0] <- sample(rep(1:nfolds, length.out = length(idx0)))
  foldid[idx1] <- sample(rep(1:nfolds, length.out = length(idx1)))
  return(foldid)
}

# ==============================================================================
# 2. SINGLE TRIAL SIMULATION WORKER
# ==============================================================================
run_single_trial_simulation <- function(sim_id, 
                                        compiled_std, 
                                        compiled_biotite, 
                                        n = 40, 
                                        skeleton = c(0.05, 0.10, 0.20, 0.35, 0.50)) {
  
  cat(sprintf("\n--- Starting Iteration %d --- \n", sim_id))
  
  # Ensure cohort has at least 6 DLT events and 6 non-events for 5-fold CV
  repeat {
    sim_data <- generate_centered_binary_scenario_data(n = n, skeleton = skeleton)
    dlt_sum  <- sum(sim_data$df_sim$dlt)
    if (dlt_sum >= 6 && dlt_sum <= (n - 6)) break
  }
  
  df_sim      <- sim_data$df_sim
  X_design    <- sim_data$X_design
  beta_true   <- sim_data$beta_true
  causal_mask <- ifelse(beta_true > 0, 1, 0)
  
  pf     <- c(rep(0, 5), rep(1, 100))
  foldid <- get_stratified_foldid(df_sim$dlt, nfolds = 5)
  
  # 1. LASSO / Elastic Net
  cat("  [1/6] Fitting LASSO & Elastic Net...\n")
  fit_lasso <- cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 1.0, penalty.factor = pf, foldid = foldid)
  fit_enet  <- cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 0.5, penalty.factor = pf, foldid = foldid)
  lasso_selected <- sum(as.matrix(coef(fit_lasso, s = "lambda.min"))[paste0("X_", 1:100), 1] != 0)
  enet_selected  <- sum(as.matrix(coef(fit_enet,  s = "lambda.min"))[paste0("X_", 1:100), 1] != 0)
  
  # 2. Ridge
  cat("  [2/6] Fitting Ridge...\n")
  fit_ridge  <- cv.glmnet(X_design, df_sim$dlt, family = "binomial", alpha = 0.0, penalty.factor = pf, foldid = foldid)
  ridge_beta <- as.matrix(coef(fit_ridge, s = "lambda.min"))[paste0("X_", 1:100), 1]
  
  # 3. Random Forest
  cat("  [3/6] Fitting Random Forest...\n")
  rf_input_df <- df_sim[, c("dose", paste0("X_", 1:100))]
  fit_rf      <- randomForest::randomForest(x = rf_input_df, y = factor(df_sim$dlt), importance = TRUE, ntree = 300)
  rf_imp_mat  <- randomForest::importance(fit_rf, type = 1)
  rf_beta     <- rf_imp_mat[paste0("X_", 1:100), ncol(rf_imp_mat)]
  
  # 4. XGBoost
  cat("  [4/6] Fitting XGBoost...\n")
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
  
  # 5. Standard Stan
  cat("  [5/6] Fitting Standard Stan...\n")
  stan_data_std <- list(Ntotal = nrow(X_design), Ncol = ncol(X_design), X = X_design, y = as.array(df_sim$dlt))
  fit_std  <- rstan::sampling(compiled_std, data = stan_data_std, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0)
  ext_std  <- rstan::extract(fit_std)
  std_beta <- apply(ext_std$betas[, 6:105], 2, mean)
  names(std_beta) <- paste0("X_", 1:100)
  
  # 6. BioTiTE V3
  cat("  [6/6] Fitting BioTiTE V3...\n")
  stan_data_bio <- list(
    Ntotal = nrow(X_design), Ncol = ncol(X_design), X = X_design, y = as.array(df_sim$dlt),
    w = rep(1.0, nrow(df_sim)), prior_means = qlogis(skeleton), prior_sd_alpha1 = 0.50, prior_sds_delta = rep(0.25, 4)
  )
  fit_bio  <- rstan::sampling(compiled_biotite, data = stan_data_bio, iter = 2000, warmup = 1000, chains = 2, cores = 1, refresh = 0)
  ext_bio  <- rstan::extract(fit_bio)
  bio_beta <- apply(ext_bio$beta_cov, 2, mean)
  names(bio_beta) <- paste0("X_", 1:100)
  
  ridge_ranks <- rank(abs(ridge_beta))
  rf_ranks    <- rank(rf_beta)
  xgb_ranks   <- rank(xgb_beta)
  std_ranks   <- rank(abs(std_beta))
  bio_ranks   <- rank(abs(bio_beta))
  
  aucs <- c(
    Ridge   = calc_feature_auc(abs(ridge_beta), causal_mask),
    RanFor  = calc_feature_auc(rf_beta, causal_mask),
    XGB     = calc_feature_auc(xgb_beta, causal_mask),
    StdStan = calc_feature_auc(abs(std_beta), causal_mask),
    BioTiTE  = calc_feature_auc(abs(bio_beta), causal_mask)
  )
  
  return(list(
    sim_id         = sim_id,
    lasso_selected = lasso_selected,
    enet_selected  = enet_selected,
    ridge_ranks    = ridge_ranks[1:5],
    rf_ranks       = rf_ranks[1:5],
    xgb_ranks      = xgb_ranks[1:5],
    std_ranks      = std_ranks[1:5],
    bio_ranks      = bio_ranks[1:5],
    aucs           = aucs
  ))
}

# ==============================================================================
# 3. SEQUENTIAL DEBUG LOOP (5 SIMULATIONS)
# ==============================================================================
cat("Compiling Stan models once on main node...\n")
compiled_std     <- rstan::stan_model(file = "binomialRegression.stan")
compiled_biotite <- rstan::stan_model(file = "BioTiTE_V3.stan")

N_debug_sims <- 5
results_list <- list()

for (s in 1:N_debug_sims) {
  results_list[[s]] <- run_single_trial_simulation(
    sim_id           = s,
    compiled_std     = compiled_std,
    compiled_biotite = compiled_biotite,
    n                = 40,
    skeleton         = c(0.05, 0.10, 0.20, 0.35, 0.50)
  )
}

# Aggregate and display
auc_mat <- t(sapply(results_list, function(x) x$aucs))
cat("\n======================================================================\n")
cat("DEBUG SUCCESSFUL: 5-SIMULATION AUC SUMMARY\n")
cat("======================================================================\n")
print(round(colMeans(auc_mat), 4))