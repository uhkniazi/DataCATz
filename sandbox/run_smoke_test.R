library(dfcrm)
library(trialr)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. WINDOWS 11 POWER LOCK (Prevent Throttling/Sleep)
# ==============================================================================
if (.Platform$OS.type == "windows") {
  cat("Locking Windows power settings to keep CPU awake during execution...\n")
  system("powercfg /change monitor-timeout-ac 0", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 0", ignore.stdout = TRUE)
}

# Clear any old progress tracking logs
if (file.exists("simulation_log.txt")) file.remove("simulation_log.txt")

# ==============================================================================
# 2. COMPILE ONCE OUTSIDE
# ==============================================================================
cat("Compiling custom Stan engine... (this takes a minute)\n")
compiled_stan_model <- rstan::stan_model(file = "TiTE_likelihood_regression.stan")

# ==============================================================================
# 3. GLOBAL EXPERIMENTAL CONFIGURATION (SMOKE TEST SETTINGS)
# ==============================================================================
set.seed(42)
N_sims     <- 6  # Testing exactly 1 batch of your 6 workers
N_patients <- 20
T_max      <- 28
target_dlt <- 0.25
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton)

# Biological Truth parameters
true_alpha <- qlogis(c(0.04, 0.10, 0.24, 0.42, 0.60)) 
true_b_sex <- 0.8
true_b_bmi <- -0.3

# Core Biostatistics Function: Calculates true individualized probability
get_true_risk <- function(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi) {
  plogis(true_alpha[dose] + (true_b_sex * female) + (true_b_bmi * bmi_centered))
}

simulate_patient_outcome <- function(dose, female, bmi_centered, arrival_day, T_max, true_alpha, true_b_sex, true_b_bmi) {
  p_true  <- get_true_risk(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  dlt     <- rbinom(1, 1, p_true)
  clearance_day <- arrival_day + ifelse(dlt == 1, sample(1:24, 1), T_max)
  return(list(dlt = dlt, clearance_day = clearance_day, p_true = p_true))
}

# ==============================================================================
# 4. PARALLEL SINGLE TRIAL FUNCTION WITH LIVE PROGRESS TIMESTAMPS
# ==============================================================================
run_single_trial_simulation <- function(sim_id, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi) {
  
  # Append progress live to the text file so you know it's not frozen
  cat(sprintf("[%s] Worker node engaged: Starting Simulation Run #%d\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  set.seed(42 + sim_id)
  
  patients <- data.frame(
    id          = 1:N_patients,
    female      = sample(c(0, 1), N_patients, replace = TRUE),
    bmi_raw     = rnorm(N_patients, mean = 26, sd = 4),
    arrival_day = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE)))
  )
  patients$bmi_centered <- patients$bmi_raw - 26
  
  d_dfcrm  <- rep(NA, N_patients); y_dfcrm  <- rep(NA, N_patients); c_dfcrm  <- rep(NA, N_patients); p_dfcrm  <- rep(NA, N_patients)
  d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients); p_trialr <- rep(NA, N_patients)
  d_m5base <- rep(NA, N_patients); y_m5base <- rep(NA, N_patients); c_m5base <- rep(NA, N_patients); p_m5base <- rep(NA, N_patients)
  d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients); p_m5cov  <- rep(NA, N_patients)
  
  for (i in 1:N_patients) {
    current_time <- patients$arrival_day[i]
    comp <- 1:(i-1)
    
    # --- DFCRM ---
    if (i == 1) { d_dfcrm[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_dfcrm[comp] == 1, c_dfcrm[comp] - patients$arrival_day[comp], days_on_study)
      w_dfcrm       <- ifelse(y_dfcrm[comp] == 1, 1.0, days_on_study / T_max)
      fit_dfcrm     <- dfcrm::titecrm(prior=q_skeleton, target=target_dlt, tox=y_dfcrm[comp], level=d_dfcrm[comp], weights=w_dfcrm, followup=T_max, model="empiric")
      d_dfcrm[i]    <- min(fit_dfcrm$mtd, d_dfcrm[i-1] + 1)
    }
    out_dfcrm <- simulate_patient_outcome(d_dfcrm[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_dfcrm[i] <- out_dfcrm$dlt; c_dfcrm[i] <- out_dfcrm$clearance_day; p_dfcrm[i] <- out_dfcrm$p_true
    
    # --- TRIALR ---
    if (i == 1) { d_trialr[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
      w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
      fit_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=1.34, chains=4, cores=4, refresh=0)
      d_trialr[i]   <- min(which.min(abs(fit_trialr$prob_tox - target_dlt)), d_trialr[i-1] + 1)
    }
    out_trialr <- simulate_patient_outcome(d_trialr[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day; p_trialr[i] <- out_trialr$p_true
    
    # --- M5 BASE ---
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
    y_m5base[i] <- out_m5base$dlt; c_m5base[i] <- out_m5base$clearance_day; p_m5base[i] <- out_m5base$p_true
    
    # --- M5 COVARIATE ---
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
    y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day; p_m5cov[i] <- out_m5cov$p_true
  }
  
  base_risks_i <- sapply(1:N_patients, function(idx) {
    get_true_risk(d_m5cov[idx], patients$female[idx], patients$bmi_centered[idx], true_alpha, true_b_sex, true_b_bmi)
  })
  
  cat(sprintf("[%s] --- Run #%d SUCCESSFULLY COMPLETED ---\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  return(list(
    mtd_dfcrm = d_dfcrm[N_patients], dlts_dfcrm = sum(y_dfcrm), patient_overdoses_dfcrm = sum(p_dfcrm > target_dlt),
    mtd_trialr = d_trialr[N_patients], dlts_trialr = sum(y_trialr), patient_overdoses_trialr = sum(p_trialr > target_dlt),
    mtd_m5base = d_m5base[N_patients], dlts_m5base = sum(y_m5base), patient_overdoses_m5base = sum(p_m5base > target_dlt),
    mtd_m5cov = d_m5cov[N_patients], dlts_m5cov = sum(y_m5cov), patient_overdoses_m5cov = sum(p_m5cov > target_dlt),
    cov_assigned_doses = d_m5cov, cov_patient_risks = base_risks_i,
    demographics = data.frame(female = patients$female, bmi = patients$bmi_raw, y_cov = y_m5cov, y_base = y_m5base)
  ))
}

# ==============================================================================
# 5. HIGH-THROUGHPUT CLUSTER PROCESSING
# ==============================================================================
num_workers <- 6
cat(sprintf("\nLaunching parallel cluster across %d workers for smoke test...\n", num_workers))
cl <- parallel::makeCluster(num_workers)
parallel::clusterEvalQ(cl, { library(dfcrm); library(trialr); library(rstan) })
parallel::clusterExport(cl, c("compiled_stan_model", "N_patients", "T_max", "target_dlt", 
                              "q_skeleton", "logit_q", "true_alpha", "true_b_sex", "true_b_bmi", 
                              "get_true_risk", "simulate_patient_outcome", "run_single_trial_simulation"))

raw_parallel_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
  run_single_trial_simulation(s, compiled_stan_model, N_patients, T_max, target_dlt, 
                              q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi)
})
parallel::stopCluster(cl)

# ==============================================================================
# 6. POST-PROCESSING SCORECARD GENERATION
# ==============================================================================
final_mtd_dfcrm  <- sapply(raw_parallel_results, function(x) x$mtd_dfcrm)
final_mtd_trialr <- sapply(raw_parallel_results, function(x) x$mtd_trialr)
final_mtd_m5base <- sapply(raw_parallel_results, function(x) x$mtd_m5base)
final_mtd_m5cov  <- sapply(raw_parallel_results, function(x) x$mtd_m5cov)

total_dlts_dfcrm  <- sapply(raw_parallel_results, function(x) x$dlts_dfcrm)
total_dlts_trialr <- sapply(raw_parallel_results, function(x) x$dlts_trialr)
total_dlts_m5base <- sapply(raw_parallel_results, function(x) x$dlts_m5base)
total_dlts_m5cov  <- sapply(raw_parallel_results, function(x) x$dlts_m5cov)

pt_overdose_dfcrm  <- sapply(raw_parallel_results, function(x) x$patient_overdoses_dfcrm)
pt_overdose_trialr <- sapply(raw_parallel_results, function(x) x$patient_overdoses_trialr)
pt_overdose_m5base <- sapply(raw_parallel_results, function(x) x$patient_overdoses_m5base)
pt_overdose_m5cov  <- sapply(raw_parallel_results, function(x) x$patient_overdoses_m5cov)

all_assigned_doses <- unlist(lapply(raw_parallel_results, function(x) x$cov_assigned_doses))
all_patient_risks   <- unlist(lapply(raw_parallel_results, function(x) x$cov_patient_risks))
risk_alignment_r    <- cor(all_patient_risks, all_assigned_doses)

all_demographics <- do.call(rbind, lapply(raw_parallel_results, function(x) x$demographics))
fragile_cohort   <- all_demographics[all_demographics$female == 1 & all_demographics$bmi < 22, ]

cat("\n========================================================================\n")
cat("          UPGRADED CLINICAL METRICS (SMOKE TEST SCORECARD)              \n")
cat("========================================================================\n")

oc_summary <- data.frame(
  Framework             = c("dfcrm (MLE)", "trialr (MCMC)", "M5Base (Pure Pop)", "M5Cov (Personalized)"),
  Global_MTD_Accuracy   = c(mean(final_mtd_dfcrm == 3), mean(final_mtd_trialr == 3), mean(final_mtd_m5base == 3), mean(final_mtd_m5cov == 3)) * 100,
  Avg_DLTs_Per_Trial    = c(mean(total_dlts_dfcrm), mean(total_dlts_trialr), mean(total_dlts_m5base), mean(total_dlts_m5cov)),
  True_Patient_Overdose = c(mean(pt_overdose_dfcrm), mean(pt_overdose_trialr), mean(pt_overdose_m5base), mean(pt_overdose_m5cov))
)
oc_summary[, 2:4] <- round(oc_summary[, 2:4], 2)
print(oc_summary)

cat("\n------------------ ADVANCED PERSONALIZATION ADVANTAGE ------------------\n")
cat(sprintf("Risk-Allocation Correlation (r) for M5Cov  : %.2f\n", risk_alignment_r))
if (nrow(fragile_cohort) > 0) {
  cat(sprintf("DLT Rate in Fragile Cohort (Low-BMI Female) under M5Base : %.1f%%\n", mean(fragile_cohort$y_base) * 100))
  cat(sprintf("DLT Rate in Fragile Cohort (Low-BMI Female) under M5Cov  : %.1f%%\n", mean(fragile_cohort$y_cov) * 100))
} else {
  cat("Fragile Cohort note: No low-BMI females happened to be sampled in this 6-run batch.\n")
}

# ==============================================================================
# 7. CLEAN UP POWER SETTINGS
# ==============================================================================
if (.Platform$OS.type == "windows") {
  cat("\nRestoring default Windows sleep timers...\n")
  system("powercfg /change monitor-timeout-ac 15", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 15", ignore.stdout = TRUE)
}