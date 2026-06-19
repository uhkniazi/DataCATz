library(dfcrm)
library(trialr)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. COMPILE ONCE OUTSIDE
# ==============================================================================
cat("Compiling your custom Stan model... (Reused natively)\n")
compiled_stan_model <- rstan::stan_model(file = "TiTE_likelihood_regression.stan")

# ==============================================================================
# 2. GLOBAL EXPERIMENTAL CONFIGURATION
# ==============================================================================
set.seed(42)
N_sims     <- 200  # Total trials to execute
N_patients <- 20
T_max      <- 28
target_dlt <- 0.25
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton)

# Biological Truth parameters
true_alpha <- qlogis(c(0.04, 0.10, 0.24, 0.42, 0.60)) # Target MTD = Dose 3
true_b_sex <- 0.8
true_b_bmi <- -0.3

# Outcome Generator Helper
simulate_patient_outcome <- function(dose, female, bmi_centered, arrival_day, T_max, true_alpha, true_b_sex, true_b_bmi) {
  logit_p <- true_alpha[dose] + (true_b_sex * female) + (true_b_bmi * bmi_centered)
  p_true  <- plogis(logit_p)
  dlt     <- rbinom(1, 1, p_true)
  clearance_day <- arrival_day + ifelse(dlt == 1, sample(1:24, 1), T_max)
  return(list(dlt = dlt, clearance_day = clearance_day))
}

# ==============================================================================
# 3. DEFINE THE INDEPENDENT SINGLE TRIAL FUNCTION
# ==============================================================================
run_single_trial_simulation <- function(sim_id, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi) {
  
  # Ensure each parallel core generates its own distinct random stream
  set.seed(42 + sim_id)
  
  patients <- data.frame(
    id          = 1:N_patients,
    female      = sample(c(0, 1), N_patients, replace = TRUE),
    bmi_raw     = rnorm(N_patients, mean = 26, sd = 4),
    arrival_day = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE)))
  )
  patients$bmi_centered <- patients$bmi_raw - 26
  
  d_dfcrm  <- rep(NA, N_patients); y_dfcrm  <- rep(NA, N_patients); c_dfcrm  <- rep(NA, N_patients)
  d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients)
  d_m5base <- rep(NA, N_patients); y_m5base <- rep(NA, N_patients); c_m5base <- rep(NA, N_patients)
  d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients)
  
  # Tracking local sub-cohort dose choices for this isolated trial
  local_m5_males    <- c(); local_m5_females  <- c()
  local_m5_low_bmi  <- c(); local_m5_high_bmi <- c()
  
  for (i in 1:N_patients) {
    current_time <- patients$arrival_day[i]
    comp <- 1:(i-1)
    
    # --- REALITY 1: DFCRM ---
    if (i == 1) { d_dfcrm[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_dfcrm[comp] == 1, c_dfcrm[comp] - patients$arrival_day[comp], days_on_study)
      w_dfcrm       <- ifelse(y_dfcrm[comp] == 1, 1.0, days_on_study / T_max)
      fit_dfcrm     <- dfcrm::titecrm(prior=q_skeleton, target=target_dlt, tox=y_dfcrm[comp], level=d_dfcrm[comp], weights=w_dfcrm, followup=T_max, model="empiric")
      d_dfcrm[i]    <- min(fit_dfcrm$mtd, d_dfcrm[i-1] + 1)
    }
    out_dfcrm <- simulate_patient_outcome(d_dfcrm[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_dfcrm[i] <- out_dfcrm$dlt; c_dfcrm[i] <- out_dfcrm$clearance_day
    
    # --- REALITY 2: TRIALR (Running 4 Chains In Parallel) ---
    if (i == 1) { d_trialr[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
      w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
      fit_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=1.34, chains=4, cores=4, refresh=0)
      d_trialr[i]   <- min(which.min(abs(fit_trialr$prob_tox - target_dlt)), d_trialr[i-1] + 1)
    }
    out_trialr <- simulate_patient_outcome(d_trialr[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day
    
    # --- REALITY 3: M5 BASE (4 Chains, 4 Cores) ---
    if (i == 1) { d_m5base[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5base[comp] == 1, c_m5base[comp] - patients$arrival_day[comp], days_on_study)
      w_m5base      <- ifelse(y_m5base[comp] == 1, 1.0, days_on_study / T_max)
      X_base        <- model.matrix(~ 0 + factor(d_m5base[comp], levels=1:5))
      stan_data_b   <- list(Ntotal=nrow(X_base), Ncol=ncol(X_base), X=X_base, y=as.array(y_m5base[comp]), w=as.array(w_m5base), prior_means=logit_q)
      fit_stan_b    <- rstan::sampling(compiled_stan_model, data=stan_data_b, iter=1000, chains=4, cores=4, refresh=0, warmup=400)
      base_alphas   <- rstan::extract(fit_stan_b)$betas
      pred_b        <- sapply(1:5, function(d) mean(plogis(base_alphas[, d])))
      d_m5base[i]   <- min(which.min(abs(pred_b - target_dlt)), d_m5base[i-1] + 1)
    }
    out_m5base <- simulate_patient_outcome(d_m5base[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5base[i] <- out_m5base$dlt; c_m5base[i] <- out_m5base$clearance_day
    
    # --- REALITY 4: M5 COVARIATE (4 Chains, 4 Cores) ---
    if (i == 1) { d_m5cov[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5cov[comp] == 1, c_m5cov[comp] - patients$arrival_day[comp], days_on_study)
      w_m5cov       <- ifelse(y_m5cov[comp] == 1, 1.0, days_on_study / T_max)
      X_cov         <- model.matrix(~ 0 + factor(d_m5cov[comp], levels=1:5) + female + bmi_centered, data=patients[comp,])
      stan_data_c   <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov[comp]), w=as.array(w_m5cov), prior_means=logit_q)
      fit_stan_c    <- rstan::sampling(compiled_stan_model, data=stan_data_c, iter=1000, chains=4, cores=4, refresh=0, warmup=400)
      cov_betas     <- rstan::extract(fit_stan_c)$betas
      pred_c        <- sapply(1:5, function(d) mean(plogis(cov_betas[, d] + (cov_betas[, 6] * patients$female[i]) + (cov_betas[, 7] * patients$bmi_centered[i]))))
      d_m5cov[i]    <- min(which.min(abs(pred_c - target_dlt)), d_m5cov[i-1] + 1)
    }
    out_m5cov <- simulate_patient_outcome(d_m5cov[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day
    
    # Collect tracking subset allocations
    if(patients$female[i] == 0) local_m5_males    <- c(local_m5_males, d_m5cov[i])
    if(patients$female[i] == 1) local_m5_females  <- c(local_m5_females, d_m5cov[i])
    if(patients$bmi_raw[i] < 22) local_m5_low_bmi  <- c(local_m5_low_bmi, d_m5cov[i])
    if(patients$bmi_raw[i] >= 28) local_m5_high_bmi <- c(local_m5_high_bmi, d_m5cov[i])
  }
  
  return(list(
    mtd_dfcrm  = d_dfcrm[N_patients],  dlts_dfcrm  = sum(y_dfcrm),  overdose_dfcrm  = sum(d_dfcrm > 3),
    mtd_trialr = d_trialr[N_patients], dlts_trialr = sum(y_trialr), overdose_trialr = sum(d_trialr > 3),
    mtd_m5base = d_m5base[N_patients], dlts_m5base = sum(y_m5base), overdose_m5base = sum(d_m5base > 3),
    mtd_m5cov  = d_m5cov[N_patients],  dlts_m5cov  = sum(y_m5cov),  overdose_m5cov  = sum(d_m5cov > 3),
    males = local_m5_males, females = local_m5_females, low_bmi = local_m5_low_bmi, high_bmi = local_m5_high_bmi
  ))
}

# ==============================================================================
# 4. PARALLEL TUNED EXECUTION (6 Workers x 4 Chains = 24 Cores Active)
# ==============================================================================
num_workers <- 6
cat(sprintf("\nSpawning %d concurrent workers on your 28-core machine...\n", num_workers))

cl <- parallel::makeCluster(num_workers)

parallel::clusterEvalQ(cl, { library(dfcrm); library(trialr); library(rstan) })
parallel::clusterExport(cl, c("compiled_stan_model", "N_patients", "T_max", "target_dlt", 
                              "q_skeleton", "logit_q", "true_alpha", "true_b_sex", "true_b_bmi", 
                              "simulate_patient_outcome", "run_single_trial_simulation"))

raw_parallel_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
  run_single_trial_simulation(s, compiled_stan_model, N_patients, T_max, target_dlt, 
                              q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi)
})

parallel::stopCluster(cl)
cat("Parallel loops successfully finished. Computing long-term statistics...\n")

# ==============================================================================
# 5. POST-PROCESSING & COMPREHENSIVE SCORECARD REPORT
# ==============================================================================
final_mtd_dfcrm  <- sapply(raw_parallel_results, function(x) x$mtd_dfcrm)
final_mtd_trialr <- sapply(raw_parallel_results, function(x) x$mtd_trialr)
final_mtd_m5base <- sapply(raw_parallel_results, function(x) x$mtd_m5base)
final_mtd_m5cov  <- sapply(raw_parallel_results, function(x) x$mtd_m5cov)

total_dlts_dfcrm  <- sapply(raw_parallel_results, function(x) x$dlts_dfcrm)
total_dlts_trialr <- sapply(raw_parallel_results, function(x) x$dlts_trialr)
total_dlts_m5base <- sapply(raw_parallel_results, function(x) x$dlts_m5base)
total_dlts_m5cov  <- sapply(raw_parallel_results, function(x) x$dlts_m5cov)

overdose_dfcrm  <- sapply(raw_parallel_results, function(x) x$overdose_dfcrm)
overdose_trialr <- sapply(raw_parallel_results, function(x) x$overdose_trialr)
overdose_m5base <- sapply(raw_parallel_results, function(x) x$overdose_m5base)
overdose_m5cov  <- sapply(raw_parallel_results, function(x) x$overdose_m5cov)

agg_males    <- unlist(lapply(raw_parallel_results, function(x) x$males))
agg_females  <- unlist(lapply(raw_parallel_results, function(x) x$females))
agg_low_bmi  <- unlist(lapply(raw_parallel_results, function(x) x$low_bmi))
agg_high_bmi <- unlist(lapply(raw_parallel_results, function(x) x$high_bmi))

cat("\n========================================================================\n")
cat("            FINAL CLINICAL TRIAL OPERATING CHARACTERISTICS              \n")
cat("========================================================================\n")

oc_summary <- data.frame(
  Framework        = c("dfcrm (MLE)", "trialr (MCMC)", "M5Base (Pure Pop)", "M5Cov (Personalized)"),
  MTD_Accuracy_Pct = c(mean(final_mtd_dfcrm == 3), mean(final_mtd_trialr == 3), mean(final_mtd_m5base == 3), mean(final_mtd_m5cov == 3)) * 100,
  Avg_DLTs_Per_Trial = c(mean(total_dlts_dfcrm), mean(total_dlts_trialr), mean(total_dlts_m5base), mean(total_dlts_m5cov)),
  Avg_Overdose_Alloc = c(mean(overdose_dfcrm), mean(overdose_trialr), mean(overdose_m5base), mean(overdose_m5cov))
)
print(round(oc_summary, 2))

cat("\n--- COVARIATE M5 PERSONALIZATION SIGNATURE CHECK ---\n")
cat(sprintf("Average Allocated Dose to Males   (Lower Risk)  : %.2f\n", mean(agg_males)))
cat(sprintf("Average Allocated Dose to Females (Higher Risk) : %.2f\n", mean(agg_females)))
cat(sprintf("Average Allocated Dose to Low BMI (<22, Fragile): %.2f\n", mean(agg_low_bmi)))
cat(sprintf("Average Allocated Dose to High BMI (>=28, Robust): %.2f\n", mean(agg_high_bmi)))