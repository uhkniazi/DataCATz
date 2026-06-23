# Name: tite_crm_simulation_A.R
# Auth: u.niazi@soton.ac.uk
# Date: 23/06/2026
# Desc: TiTE-CRM Models used in Phase 1 clinical trials to assess dosage based
#       on toxicity skeleton. Simulation for Scenario A (The "Do No Harm" Check: 
#       Covariates are completely irrelevant, skeleton is globally correct)

library(dfcrm)
library(trialr)
library(rstan)
library(parallel)

# Optimize Stan execution by writing the compiled C++ model to disk to avoid recompiling
rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. WINDOWS 11 POWER LOCK (Prevent CPU Throttling or System Sleep)
# ==============================================================================
# [Logic Decision] Parallel clusters thrashing multiple CPU cores can cause 
# Windows to aggressively throttle or trigger sleep modes during long executions.
if (.Platform$OS.type == "windows") {
  cat("Locking Windows power settings to keep CPU awake during execution...\n")
  # Set monitor and system sleep timeouts to 0 (Never Sleep) while running
  system("powercfg /change monitor-timeout-ac 0", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 0", ignore.stdout = TRUE)
}

# Clean slate: Erase old multi-core execution logs from previous sessions
if (file.exists("simulation_log.txt")) file.remove("simulation_log.txt")

# ==============================================================================
# 2. COMPILE STAN ENGINE (Executed once outside the parallel loop)
# ==============================================================================
# [Logic Decision] Compiling a Stan model takes 1-2 minutes. We compile it once 
# here on the master process and pass the compiled object to the worker nodes, 
# preventing the workers from attempting to compile it simultaneously.
cat("Compiling custom Monotonic Stan engine... (this takes a minute)\n")
compiled_stan_model <- rstan::stan_model(file = "TiTE_regression_monotonic_skeleton.stan")

# ==============================================================================
# 3. GLOBAL EXPERIMENTAL CONFIGURATION
# ==============================================================================
set.seed(42)
N_sims     <- 200  # Rigorous evaluation size to ensure statistical stability
N_patients <- 20   # Sample size per individual clinical trial
T_max      <- 28   # DLT evaluation window (28 days)
target_dlt <- 0.25 # The Maximum Tolerated Dose target (25% toxicity rate)

# The Clinical Skeleton: Standard initial population-level risk assumptions
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton) # Convert skeleton risks to logit scale for Stan priors

# --- BIOLOGICAL TRUTH PARAMETERS (Scenario A: Covariates have ZERO effect) ---
true_alpha <- qlogis(c(0.05, 0.12, 0.25, 0.40, 0.55)) 
true_b_sex <- 0.0   # [Scenario A Baseline] Females have 0.0 log-odds change (No Effect)
true_b_bmi <- 0.0   # [Scenario A Baseline] BMI changes shift log-odds by 0.0 (No Effect)

# [Mathematical Formula] Calculate a patient's true biological probability of a DLT
get_true_risk <- function(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi) {
  plogis(true_alpha[dose] + (true_b_sex * female) + (true_b_bmi * bmi_centered))
}

# [Clinical Decision] Identify the exact dose closest to the 25% target for a specific profile
get_true_optimal_dose <- function(female, bmi_centered, true_alpha, true_b_sex, true_b_bmi, target_dlt = 0.25) {
  true_risks <- sapply(1:5, function(d) {
    get_true_risk(d, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  })
  # Returns the dose level (1-5) that minimizes absolute distance to target_dlt
  return(which.min(abs(true_risks - target_dlt)))
}

# [Trial Simulator] Simulates real-time patient outcomes and delayed DLT time-to-events
simulate_patient_outcome <- function(dose, female, bmi_centered, arrival_day, T_max, true_alpha, true_b_sex, true_b_bmi) {
  p_true  <- get_true_risk(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  dlt     <- rbinom(1, 1, p_true) # Coin flip using the patient's true uniform population probability
  
  # If a DLT occurs, it happens randomly between days 1 and 24. If not, they clear T_max (28 days)
  clearance_day <- arrival_day + ifelse(dlt == 1, sample(1:24, 1), T_max)
  return(list(dlt = dlt, clearance_day = clearance_day, p_true = p_true))
}

# ==============================================================================
# 4. PARALLEL SINGLE TRIAL FUNCTION (Executed simultaneously across workers)
# ==============================================================================
run_single_trial_simulation <- function(sim_id, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi) {
  
  # Log progress to disk so the user can monitor high-throughput execution in real time
  cat(sprintf("[%s] Worker node engaged: Starting Simulation Run #%d\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  # Ensure strict multi-core reproducibility by isolating the random seed on each worker
  set.seed(42 + sim_id)
  
  # Generate heterogeneous patient demographics sequentially entering the clinic
  patients <- data.frame(
    id          = 1:N_patients,
    female      = sample(c(0, 1), N_patients, replace = TRUE),
    bmi_raw     = rnorm(N_patients, mean = 26, sd = 4),
    arrival_day = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE))) # Staggered sequential accrual
  )
  # Center BMI around the population average (26) to make regression intercepts interpretable
  patients$bmi_centered <- patients$bmi_raw - 26
  
  # Calculate the exact personalized dose target for every patient entering this specific trial
  # In Scenario A, because beta=0, every single patient's true optimal dose will safely collapse to Dose 3.
  patients$true_optimal <- sapply(1:N_patients, function(idx) {
    get_true_optimal_dose(patients$female[idx], patients$bmi_centered[idx], true_alpha, true_b_sex, true_b_bmi, target_dlt)
  })
  
  # Allocate tracking arrays for the four competing statistical frameworks
  d_dfcrm  <- rep(NA, N_patients); y_dfcrm  <- rep(NA, N_patients); c_dfcrm  <- rep(NA, N_patients); p_dfcrm  <- rep(NA, N_patients)
  d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients); p_trialr <- rep(NA, N_patients)
  d_m5base <- rep(NA, N_patients); y_m5base <- rep(NA, N_patients); c_m5base <- rep(NA, N_patients); p_m5base <- rep(NA, N_patients)
  d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients); p_m5cov  <- rep(NA, N_patients)
  
  # --- BEGIN SEQUENTIAL PATIENT ACCRUAL LOOP ---
  for (i in 1:N_patients) {
    current_time <- patients$arrival_day[i] # Current calendar day of the trial
    comp <- 1:(i-1)                         # Index of all historical patients enrolled so far
    
    # ==========================================================================
    # FRAMEWORK 1: CLASSIC MAXIMUM LIKELIHOOD ESTIMATION (dfcrm)
    # ==========================================================================
    # [Clinical Decision] The first patient in a Phase I trial must always be assigned 
    # to Dose Level 1 for safety. We lack data to authorize a higher step.
    if (i == 1) { d_dfcrm[i] <- 1 } else {
      # [TiTE Calculus] Compute time-on-study for previous patients relative to current day
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      # If a patient had a DLT, lock their follow-up duration to their exact toxicity event day
      days_on_study <- ifelse(y_dfcrm[comp] == 1, c_dfcrm[comp] - patients$arrival_day[comp], days_on_study)
      # Calculate partial follow-up weights: Complete follow-up or DLT = 1.0; Partial follow-up = fraction of T_max
      w_dfcrm       <- ifelse(y_dfcrm[comp] == 1, 1.0, days_on_study / T_max)
      
      fit_dfcrm     <- dfcrm::titecrm(prior=q_skeleton, target=target_dlt, tox=y_dfcrm[comp], level=d_dfcrm[comp], weights=w_dfcrm, followup=T_max, model="empiric")
      # [Clinical Safety Rule] Force the "No-Skipping" constraint: The model can never 
      # escalate higher than 1 dose level above the last treated patient's dose.
      d_dfcrm[i]    <- min(fit_dfcrm$mtd, d_dfcrm[i-1] + 1)
    }
    out_dfcrm <- simulate_patient_outcome(d_dfcrm[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_dfcrm[i] <- out_dfcrm$dlt; c_dfcrm[i] <- out_dfcrm$clearance_day; p_dfcrm[i] <- out_dfcrm$p_true
    
    # ==========================================================================
    # FRAMEWORK 2: BAYESIAN SINGLE-PARAMETER MCMC (trialr)
    # ==========================================================================
    # [Clinical Decision] Enforce safety baseline: Patient 1 starts at Dose 1.
    if (i == 1) { d_trialr[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
      w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Prior Adjustment] Pass beta_sd = sqrt(1.34) to match standard historical prior variance of 1.34
      fit_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=sqrt(1.34), chains=4, cores=1, refresh=0)
      # [Clinical Safety Rule] Enforce the "No-Skipping" escalation cap
      d_trialr[i]   <- min(which.min(abs(fit_trialr$prob_tox - target_dlt)), d_trialr[i-1] + 1)
    }
    out_trialr <- simulate_patient_outcome(d_trialr[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day; p_trialr[i] <- out_trialr$p_true
    
    # ==========================================================================
    # FRAMEWORK 3: CUSTOM PURE POPULATION MODEL (M5Base)
    # ==========================================================================
    # [Clinical Decision] Enforce safety baseline: Patient 1 starts at Dose 1.
    if (i == 1) { d_m5base[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5base[comp] == 1, c_m5base[comp] - patients$arrival_day[comp], days_on_study)
      w_m5base      <- ifelse(y_m5base[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Pseudocode] Build structural design matrix using strictly population dose indicators (blind to patient traits)
      X_base        <- model.matrix(~ 0 + factor(d_m5base[comp], levels=1:5))
      stan_data_b   <- list(Ntotal=nrow(X_base), Ncol=ncol(X_base), X=X_base, y=as.array(y_m5base[comp]), w=as.array(w_m5base), prior_means=logit_q)
      fit_stan_b    <- rstan::sampling(compiled_stan_model, data=stan_data_b, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
      
      base_alphas   <- rstan::extract(fit_stan_b)$alpha  
      pred_b        <- sapply(1:5, function(d) mean(plogis(base_alphas[, d]))) # Mean posterior toxicity curve
      # [Clinical Safety Rule] Enforce the "No-Skipping" escalation cap
      d_m5base[i]   <- min(which.min(abs(pred_b - target_dlt)), d_m5base[i-1] + 1)
    }
    out_m5base <- simulate_patient_outcome(d_m5base[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5base[i] <- out_m5base$dlt; c_m5base[i] <- out_m5base$clearance_day; p_m5base[i] <- out_m5base$p_true
    
    # ==========================================================================
    # FRAMEWORK 4: CUSTOM COVARIATE-PERSONALIZED MODEL (M5Cov)
    # ==============================================================================
    # [Clinical Decision] Enforce safety baseline: Patient 1 starts at Dose 1.
    if (i == 1) { d_m5cov[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5cov[comp] == 1, c_m5cov[comp] - patients$arrival_day[comp], days_on_study)
      w_m5cov       <- ifelse(y_m5cov[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Pseudocode] Build complete design matrix tracking both Dose Levels AND individual Patient Covariates
      X_cov         <- model.matrix(~ 0 + factor(d_m5cov[comp], levels=1:5) + female + bmi_centered, data=patients[comp,])
      stan_data_c   <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov[comp]), w=as.array(w_m5cov), prior_means=logit_q)
      fit_stan_c    <- rstan::sampling(compiled_stan_model, data=stan_data_c, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
      
      extracted_c   <- rstan::extract(fit_stan_c)
      cov_alpha     <- extracted_c$alpha     
      cov_beta_cov  <- extracted_c$beta_cov  
      
      # [Personalization Logic] Project the posterior toxicity curve tailored to the i-th patient's profile.
      # In Scenario A, the model should ideally learn that cov_beta_cov values hover near 0.
      pred_c        <- sapply(1:5, function(d) {
        mean(plogis(cov_alpha[, d] + (cov_beta_cov[, 1] * patients$female[i]) + (cov_beta_cov[, 2] * patients$bmi_centered[i])))
      })
      # [Clinical Safety Rule] Enforce the "No-Skipping" escalation cap
      d_m5cov[i]    <- min(which.min(abs(pred_c - target_dlt)), d_m5cov[i-1] + 1)
    }
    out_m5cov <- simulate_patient_outcome(d_m5cov[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day; p_m5cov[i] <- out_m5cov$p_true
  }
  
  # Complete execution verification statement for multi-core tracking logs
  cat(sprintf("[%s] --- Run #%d SUCCESSFULLY COMPLETED ---\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  # Return data structures containing metrics computed across the simulated patient sample
  return(list(
    # [Final MTD Selections] The final dose assignment recommended for the *last* patient (Patient 20)
    mtd_dfcrm = d_dfcrm[N_patients], mtd_trialr = d_trialr[N_patients], mtd_m5base = d_m5base[N_patients], mtd_m5cov = d_m5cov[N_patients],
    
    # [Toxicity Burden] Total observed DLT events across the 20 patients in this trial
    dlts_dfcrm = sum(y_dfcrm), dlts_trialr = sum(y_trialr), dlts_m5base = sum(y_m5base), dlts_m5cov = sum(y_m5cov),
    
    # [Overdose Metrics] Counts how many individual patients were given a dose where their 
    # true biological risk exceeded the target (risk > 0.25). Essential for safety auditing.
    patient_overdoses_dfcrm  = sum(p_dfcrm > target_dlt),  patient_overdoses_trialr = sum(p_trialr > target_dlt),
    patient_overdoses_m5base = sum(p_m5base > target_dlt), patient_overdoses_m5cov  = sum(p_m5cov > target_dlt),
    
    # [Exact Personalization Metric] The percentage of patients who received their exact 
    # personalized optimal dose level. Shows perfect-match precision.
    personal_acc_dfcrm  = mean(d_dfcrm == patients$true_optimal) * 100,
    personal_acc_trialr = mean(d_trialr == patients$true_optimal) * 100,
    personal_acc_m5base = mean(d_m5base == patients$true_optimal) * 100,
    personal_acc_m5cov  = mean(d_m5cov == patients$true_optimal) * 100,
    
    # [Clinical Margin Metric] Plus-or-Minus 1 Dose Level Accuracy. Calculates how often 
    # a model assigned a dose safely within one level of the optimal point.
    personal_pm1_dfcrm  = mean(abs(d_dfcrm - patients$true_optimal) <= 1) * 100,
    personal_pm1_trialr = mean(abs(d_trialr - patients$true_optimal) <= 1) * 100,
    personal_pm1_m5base = mean(abs(d_m5base - patients$true_optimal) <= 1) * 100,
    personal_pm1_m5cov  = mean(abs(d_m5cov - patients$true_optimal) <= 1) * 100,
    
    # Pool demographic profiles to isolate the vulnerable patient subgroup (Low-BMI Females)
    demographics = data.frame(female = patients$female, bmi = patients$bmi_raw, y_cov = y_m5cov, y_base = y_m5base)
  ))
}

# ==============================================================================
# 5. HIGH-THROUGHPUT PARALLEL CLUSTER PROCESSING
# ==============================================================================
# Set up worker count: Consume maximum available processing threads minus 4 to prevent system lockups
num_workers <- min(detectCores() - 4, N_sims)
cat(sprintf("\nLaunching parallel cluster across %d workers for Scenario A...\n", num_workers))

cl <- parallel::makeCluster(num_workers)
# Load libraries inside each separate sub-process environment
parallel::clusterEvalQ(cl, { library(dfcrm); library(trialr); library(rstan) })
# Export configurations, definitions, and functions to the parallel workspace
parallel::clusterExport(cl, c("compiled_stan_model", "N_patients", "T_max", "target_dlt", 
                              "q_skeleton", "logit_q", "true_alpha", "true_b_sex", "true_b_bmi", 
                              "get_true_risk", "get_true_optimal_dose", "simulate_patient_outcome", "run_single_trial_simulation"))

# [Pseudocode] Execute the loop simultaneously across workers and collect results
raw_parallel_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
  run_single_trial_simulation(s, compiled_stan_model, N_patients, T_max, target_dlt, 
                              q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi)
})
parallel::stopCluster(cl) # Clean up and close the cluster connections

# ==============================================================================
# 6. POST-PROCESSING SUMMARY SCORECARD GENERATION
# ==============================================================================
# Extract metrics from all across-trial lists for structural summary calculation
final_mtd_dfcrm  <- sapply(raw_parallel_results, function(x) x$mtd_dfcrm)
final_mtd_trialr <- sapply(raw_parallel_results, function(x) x$mtd_trialr)
final_mtd_m5base <- sapply(raw_parallel_results, function(x) x$mtd_m5base)

total_dlts_dfcrm  <- sapply(raw_parallel_results, function(x) x$dlts_dfcrm)
total_dlts_trialr <- sapply(raw_parallel_results, function(x) x$dlts_trialr)
total_dlts_m5base <- sapply(raw_parallel_results, function(x) x$dlts_m5base)
total_dlts_m5cov  <- sapply(raw_parallel_results, function(x) x$dlts_m5cov)

pt_overdose_dfcrm  <- sapply(raw_parallel_results, function(x) x$patient_overdoses_dfcrm)
pt_overdose_trialr <- sapply(raw_parallel_results, function(x) x$patient_overdoses_trialr)
pt_overdose_m5base <- sapply(raw_parallel_results, function(x) x$patient_overdoses_m5base)
pt_overdose_m5cov  <- sapply(raw_parallel_results, function(x) x$patient_overdoses_m5cov)

p_acc_dfcrm  <- sapply(raw_parallel_results, function(x) x$personal_acc_dfcrm)
p_acc_trialr <- sapply(raw_parallel_results, function(x) x$personal_acc_trialr)
p_acc_m5base <- sapply(raw_parallel_results, function(x) x$personal_acc_m5base)
p_acc_m5cov  <- sapply(raw_parallel_results, function(x) x$personal_acc_m5cov)

p_pm1_dfcrm  <- sapply(raw_parallel_results, function(x) x$personal_pm1_dfcrm)
p_pm1_trialr <- sapply(raw_parallel_results, function(x) x$personal_pm1_trialr)
p_pm1_m5base <- sapply(raw_parallel_results, function(x) x$personal_pm1_m5base)
p_pm1_m5cov  <- sapply(raw_parallel_results, function(x) x$personal_pm1_m5cov)

# Flatten and bind demographic tables to evaluate high-risk patient subgroups
all_demographics <- do.call(rbind, lapply(raw_parallel_results, function(x) x$demographics))
# [Subgroup Isolation] Filter patients matching high-risk criteria: Female AND Low BMI (< 22)
fragile_cohort   <- all_demographics[all_demographics$female == 1 & all_demographics$bmi < 22, ]

cat("\n========================================================================\n")
cat("          PRODUCTION SCORECARD: SCENARIO A EVALUATION (DO NO HARM)      \n")
cat("========================================================================\n")

# [Pseudocode] Construct summary data frame, average figures across runs, and print output
oc_summary <- data.frame(
  Framework             = c("dfcrm (MLE)", "trialr (MCMC)", "M5Base (Pure Pop)", "M5Cov (Personalized)"),
  
  # Global MTD Accuracy measures how often the last assigned dose is exactly Dose 3.
  # Note: Left as NA for M5Cov because a personalized model is conceptually designed to 
  # bypass population-level metrics.
  Global_MTD_Accuracy   = c(mean(final_mtd_dfcrm == 3), mean(final_mtd_trialr == 3), mean(final_mtd_m5base == 3), NA) * 100,
  
  Exact_Personal_Acc    = c(mean(p_acc_dfcrm), mean(p_acc_trialr), mean(p_acc_m5base), mean(p_acc_m5cov)),
  PlusMinus1_Dose_Acc   = c(mean(p_pm1_dfcrm), mean(p_pm1_trialr), mean(p_pm1_m5base), mean(p_pm1_m5cov)),
  Avg_DLTs_Per_Trial    = c(mean(total_dlts_dfcrm), mean(total_dlts_trialr), mean(total_dlts_m5base), mean(total_dlts_m5cov)),
  Avg_Patient_Overdoses = c(mean(pt_overdose_dfcrm), mean(pt_overdose_trialr), mean(pt_overdose_m5base), mean(pt_overdose_m5cov))
)
# Round metrics to 2 decimal places for clean, presentation-ready formatting
oc_summary[, 2:6] <- round(oc_summary[, 2:6], 2)
print(oc_summary)

cat("\n------------------ SCENARIO A PROOF OF EQUIVALENCE ---------------------\n")
if (nrow(fragile_cohort) > 0) {
  # In Scenario A, these rates should converge closely, proving the covariate adjustment safely "turns off"
  cat(sprintf("DLT Rate in Fragile Cohort (Low-BMI Female) under M5Base : %.1f%%\n", mean(fragile_cohort$y_base) * 100))
  cat(sprintf("DLT Rate in Fragile Cohort (Low-BMI Female) under M5Cov  : %.1f%%\n", mean(fragile_cohort$y_cov) * 100))
}

# ==============================================================================
# 7. CLEAN UP POWER SETTINGS
# ==============================================================================
if (.Platform$OS.type == "windows") {
  # Restore default power management schemes so your computer sleeps normally when left idle
  cat("\nRestoring default Windows sleep timers...\n")
  system("powercfg /change monitor-timeout-ac 15", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 15", ignore.stdout = TRUE)
}
