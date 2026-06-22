# ==============================================================================
# LIBRARIES & INITIALIZATION
# ==============================================================================
# [Human] Load packages for standard clinical trial design, Bayesian MCMC modeling, and parallel CPU computing.
library(dfcrm)   # Classic maximum likelihood Time-to-Event CRM
library(trialr)  # Standard Bayesian CRM package (for comparison)
library(rstan)   # Interface to Stan MCMC sampler for our custom models
library(parallel)# Spawns background CPU processes to run simulations simultaneously

# [Pseudocode] Avoid rebuilding the model binaries if the source code hasn't changed.
# [Human] Prevents R from spending a minute recompiling the Stan code every single time you hit execute.
rstan_options(auto_write = TRUE)


# ==============================================================================
# 1. OS-LEVEL POWER OPTIMIZATION
# ==============================================================================
if (.Platform$OS.type == "windows") {
  # [Pseudocode] Call system shell execution to modify OS power profiles.
  # [Human] Keeps Windows from putting a CPU core to sleep or down-throttling during long multi-hour simulations.
  cat("Locking Windows power settings to keep CPU awake during execution...\n")
  system("powercfg /change monitor-timeout-ac 0", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 0", ignore.stdout = TRUE)
}

# [Pseudocode] Delete file if exists("simulation_log.txt")
# [Human] Wipes out old logs so you have a fresh progress sheet for this run.
if (file.exists("simulation_log.txt")) file.remove("simulation_log.txt")


# ==============================================================================
# 2. STAN COMPILATION (COMPILE ONCE)
# ==============================================================================
# [Pseudocode] Parse, check syntax, translate to C++, and compile our Stan code into a machine-code object.
# [Human] We compile the model out here ONCE. This compiled engine is then passed intact to all your cluster workers.
cat("Compiling custom Monotonic Stan engine... (this takes a minute)\n")
compiled_stan_model <- rstan::stan_model(file = "TiTE_regression_monotonic_skeleton.stan")


# ==============================================================================
# 3. GLOBAL EXPERIMENTAL CONFIGURATION (SMOKE TEST SETTINGS)
# ==============================================================================
# [Pseudocode] Seed the random number generator for global reproducibility.
set.seed(42)

# [Human] Structural simulation bounds. We simulate 6 trials total, with 20 patients each.
N_sims     <- 24   # Run 6 trials (perfectly matches 1 batch on a 6-core processor)
N_patients <- 20  # The sequential sample size of each simulated trial
T_max      <- 28  # Follow-up evaluation window per patient (e.g., 28 days)
target_dlt <- 0.25# Target Dose-Limiting Toxicity rate (we want to find doses closest to 25% risk)

# [Pseudocode] q_skeleton = vector of clinical guesses. logit_q = log(q / (1 - q))
# [Human] The "Skeleton" is the clinicians' baseline guess of toxicity risk across the 5 dose levels.
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton) # Convert probabilities to logit scale for Stan priors

# [Pseudocode] Define "Ground Truth" biological simulation parameters.
# [Human] This is the underlying reality hidden from the models. The simulation engine uses this to generate data.
# Scenario B: True baseline intercepts perfectly match our skeleton.
true_alpha <- qlogis(c(0.05, 0.12, 0.25, 0.40, 0.55)) 
true_b_sex <- 0.8   # Being female adds +0.8 to logit toxicity risk (increases risk)
true_b_bmi <- -0.3  # Every unit increase in centered BMI subtracts -0.3 risk (protective effect)

# [Pseudocode] Return: 1 / (1 + exp(- (alpha[dose] + b1*Sex + b2*BMI)))
# [Human] Computes the absolute, true biological probability of toxicity for a specific patient at a specific dose.
get_true_risk <- function(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi) {
  plogis(true_alpha[dose] + (true_b_sex * female) + (true_b_bmi * bmi_centered))
}

# [Pseudocode] Scan doses 1-5, calculate true risk for this patient, return index of minimum absolute distance to 0.25
# [Human] Calculates the exact, personalized "Maximum Tolerated Dose" for an individual patient based on their traits.
get_true_optimal_dose <- function(female, bmi_centered, true_alpha, true_b_sex, true_b_bmi, target_dlt = 0.25) {
  true_risks <- sapply(1:5, function(d) {
    get_true_risk(d, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  })
  return(which.min(abs(true_risks - target_dlt))) # Returns dose number (1 to 5)
}

# [Pseudocode] Flip a biased coin based on true risk to see if DLT occurs. If yes, pick a random day for the event.
# [Human] Simulates the actual clinical outcome of a patient. If they toxicity, they trigger a DLT on a random day.
simulate_patient_outcome <- function(dose, female, bmi_centered, arrival_day, T_max, true_alpha, true_b_sex, true_b_bmi) {
  p_true  <- get_true_risk(dose, female, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  dlt     <- rbinom(1, 1, p_true) # Biased coin flip (1 = Toxicity event, 0 = Clean)
  
  # If DLT occurs, it happens early (days 1-24). If no DLT, study follow-up clears out at full T_max.
  clearance_day <- arrival_day + ifelse(dlt == 1, sample(1:24, 1), T_max)
  return(list(dlt = dlt, clearance_day = clearance_day, p_true = p_true))
}


# ==============================================================================
# 4. THE CORE SIMULATION LOOP (RUNS ONE COMPLETE TRIAL)
# ==============================================================================
run_single_trial_simulation <- function(sim_id, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi) {
  
  # [Human] Write a live timestamp to disk so you can watch execution progress via the log file.
  cat(sprintf("[%s] Worker node engaged: Starting Simulation Run #%d\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  # [Pseudocode] Establish isolated local seed per parallel worker thread.
  set.seed(42 + sim_id)
  
  # [Human] Generate 20 patients arriving sequentially with varied background characteristics (Sex, BMI).
  patients <- data.frame(
    id          = 1:N_patients,
    female      = sample(c(0, 1), N_patients, replace = TRUE),                     # 50/50 sex split
    bmi_raw     = rnorm(N_patients, mean = 26, sd = 4),                            # Standard BMI distribution
    arrival_day = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE)))       # Accrual spacing (new patient every 4-10 days)
  )
  patients$bmi_centered <- patients$bmi_raw - 26 # Center BMI at the population average of 26
  
  # [Pseudocode] Pre-calculate ground-truth benchmarks for personalization.
  patients$true_optimal <- sapply(1:N_patients, function(idx) {
    get_true_optimal_dose(patients$female[idx], patients$bmi_centered[idx], true_alpha, true_b_sex, true_b_bmi, target_dlt)
  })
  
  # [Human] Instantiate blank vectors to log assigned doses, DLTs, clearance schedules, and true risks across models.
  d_dfcrm  <- rep(NA, N_patients); y_dfcrm  <- rep(NA, N_patients); c_dfcrm  <- rep(NA, N_patients); p_dfcrm  <- rep(NA, N_patients)
  d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients); p_trialr <- rep(NA, N_patients)
  d_m5base <- rep(NA, N_patients); y_m5base <- rep(NA, N_patients); c_m5base <- rep(NA, N_patients); p_m5base <- rep(NA, N_patients)
  d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients); p_m5cov  <- rep(NA, N_patients)
  
  # ----------------------------------------------------------------------------
  # SEQUENTIAL PATIENT ACCRUAL LOOP
  # ----------------------------------------------------------------------------
  for (i in 1:N_patients) {
    current_time <- patients$arrival_day[i] # Current trial clock time when patient 'i' steps into the clinic
    comp <- 1:(i-1)                         # Index vector pointing to all previous enrolled patients
    
    # ==========================================
    # --- MODEL 1: DFCRM (MLE TiTE-CRM) --------
    # ==========================================
    if (i == 1) { d_dfcrm[i] <- 1 } else {
      # Pseudocode: w = min(current_time - arrival, T_max) / T_max. If DLT, w = 1.0
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_dfcrm[comp] == 1, c_dfcrm[comp] - patients$arrival_day[comp], days_on_study)
      w_dfcrm       <- ifelse(y_dfcrm[comp] == 1, 1.0, days_on_study / T_max) # TiTE linear fractional weights
      
      # Execute maximum likelihood estimation
      fit_dfcrm     <- dfcrm::titecrm(prior=q_skeleton, target=target_dlt, tox=y_dfcrm[comp], level=d_dfcrm[comp], weights=w_dfcrm, followup=T_max, model="empiric")
      d_dfcrm[i]    <- min(fit_dfcrm$mtd, d_dfcrm[i-1] + 1) # Assign calculated MTD, capping escalation jump at +1 dose level
    }
    # Simulate outcome for this patient under DFCRM selection
    out_dfcrm <- simulate_patient_outcome(d_dfcrm[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_dfcrm[i] <- out_dfcrm$dlt; c_dfcrm[i] <- out_dfcrm$clearance_day; p_dfcrm[i] <- out_dfcrm$p_true
    
    # ==========================================
    # --- MODEL 2: TRIALR (Bayesian TiTE-CRM) --
    # ==========================================
    if (i == 1) { d_trialr[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
      w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Human] Cores = 1 is locked here. The trialr package runs standard Bayesian CRM with an empirical power model.
      fit_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=sqrt(1.34), chains=4, cores=1, refresh=0)
      d_trialr[i]   <- min(which.min(abs(fit_trialr$prob_tox - target_dlt)), d_trialr[i-1] + 1)
    }
    out_trialr <- simulate_patient_outcome(d_trialr[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day; p_trialr[i] <- out_trialr$p_true
    
    # ==========================================
    # --- MODEL 3: M5BASE (Monotonic Pop Model)
    # ==========================================
    if (i == 1) { d_m5base[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5base[comp] == 1, c_m5base[comp] - patients$arrival_day[comp], days_on_study)
      w_m5base      <- ifelse(y_m5base[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Pseudocode] Build a clean 5-column indicator design matrix based on assigned doses.
      X_base        <- model.matrix(~ 0 + factor(d_m5base[comp], levels=1:5))
      stan_data_b   <- list(Ntotal=nrow(X_base), Ncol=ncol(X_base), X=X_base, y=as.array(y_m5base[comp]), w=as.array(w_m5base), prior_means=logit_q)
      
      # [Human] Run our custom Stan model without covariates. (Cores = 1 keeps it fast inside the cluster worker thread)
      fit_stan_b    <- rstan::sampling(compiled_stan_model, data=stan_data_b, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
      
      # [Pseudocode] Extract absolute intercepts vector. Map posterior mean through inverse-logit transformation.
      base_alphas   <- rstan::extract(fit_stan_b)$alpha  
      pred_b        <- sapply(1:5, function(d) mean(plogis(base_alphas[, d])))
      d_m5base[i]   <- min(which.min(abs(pred_b - target_dlt)), d_m5base[i-1] + 1)
    }
    out_m5base <- simulate_patient_outcome(d_m5base[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5base[i] <- out_m5base$dlt; c_m5base[i] <- out_m5base$clearance_day; p_m5base[i] <- out_m5base$p_true
    
    # ==========================================
    # --- MODEL 4: M5COV (Covariate-Adjusted) --
    # ==========================================
    if (i == 1) { d_m5cov[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5cov[comp] == 1, c_m5cov[comp] - patients$arrival_day[comp], days_on_study)
      w_m5cov       <- ifelse(y_m5cov[comp] == 1, 1.0, days_on_study / T_max)
      
      # [Pseudocode] Expand design matrix to 7 columns: [5 Dose Factors] + [Sex] + [Centered BMI]
      X_cov         <- model.matrix(~ 0 + factor(d_m5cov[comp], levels=1:5) + female + bmi_centered, data=patients[comp,])
      stan_data_c   <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov[comp]), w=as.array(w_m5cov), prior_means=logit_q)
      
      # Run MCMC sampling to generate posterior distributions for the intercepts and regression coefficients
      fit_stan_c    <- rstan::sampling(compiled_stan_model, data=stan_data_c, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
      
      extracted_c   <- rstan::extract(fit_stan_c)
      cov_alpha     <- extracted_c$alpha     # Matrix of posterior samples for the dose baseline curves
      cov_beta_cov  <- extracted_c$beta_cov  # Matrix of posterior samples for patient clinical covariates
      
      # [Human] Compute individualized toxicity predictions *specifically* for patient 'i' based on their attributes.
      pred_c        <- sapply(1:5, function(d) {
        mean(plogis(cov_alpha[, d] + (cov_beta_cov[, 1] * patients$female[i]) + (cov_beta_cov[, 2] * patients$bmi_centered[i])))
      })
      d_m5cov[i]    <- min(which.min(abs(pred_c - target_dlt)), d_m5cov[i-1] + 1)
    }
    out_m5cov <- simulate_patient_outcome(d_m5cov[i], patients$female[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day; p_m5cov[i] <- out_m5cov$p_true
  }
  
  # ----------------------------------------------------------------------------
  # TRIAL DATA WRAP-UP & PACKAGING
  # ----------------------------------------------------------------------------
  base_risks_i <- sapply(1:N_patients, function(idx) {
    get_true_risk(d_m5cov[idx], patients$female[idx], patients$bmi_centered[idx], true_alpha, true_b_sex, true_b_bmi)
  })
  
  cat(sprintf("[%s] --- Run #%d SUCCESSFULLY COMPLETED ---\n", Sys.time(), sim_id), 
      file = "simulation_log.txt", append = TRUE)
  
  # [Pseudocode] Return structured list of summary metrics for this specific trial run to the cluster collector.
  return(list(
    mtd_dfcrm = d_dfcrm[N_patients], dlts_dfcrm = sum(y_dfcrm), patient_overdoses_dfcrm = sum(p_dfcrm > target_dlt),
    mtd_trialr = d_trialr[N_patients], dlts_trialr = sum(y_trialr), patient_overdoses_trialr = sum(p_trialr > target_dlt),
    mtd_m5base = d_m5base[N_patients], dlts_m5base = sum(y_m5base), patient_overdoses_m5base = sum(p_m5base > target_dlt),
    mtd_m5cov = d_m5cov[N_patients], dlts_m5cov = sum(y_m5cov), patient_overdoses_m5cov = sum(p_m5cov > target_dlt),
    
    # Calculate Personalization Accuracy: What % of individual allocations perfectly matched their true biological optimal dose?
    personal_acc_dfcrm  = mean(d_dfcrm == patients$true_optimal) * 100,
    personal_acc_trialr = mean(d_trialr == patients$true_optimal) * 100,
    personal_acc_m5base = mean(d_m5base == patients$true_optimal) * 100,
    personal_acc_m5cov  = mean(d_m5cov == patients$true_optimal) * 100,
    
    cov_assigned_doses = d_m5cov, cov_patient_risks = base_risks_i,
    demographics = data.frame(female = patients$female, bmi = patients$bmi_raw, y_cov = y_m5cov, y_base = y_m5base)
  ))
}


# ==============================================================================
# 5. HIGH-THROUGHPUT CLUSTER PROCESSING
# ==============================================================================
num_workers <- 24
cat(sprintf("\nLaunching parallel cluster across %d workers for smoke test...\n", num_workers))

# [Pseudocode] cl = Spawn background vanilla R sessions. Export variables/compiled binaries to background memory.
cl <- parallel::makeCluster(num_workers)
parallel::clusterEvalQ(cl, { library(dfcrm); library(trialr); library(rstan) })
parallel::clusterExport(cl, c("compiled_stan_model", "N_patients", "T_max", "target_dlt", 
                              "q_skeleton", "logit_q", "true_alpha", "true_b_sex", "true_b_bmi", 
                              "get_true_risk", "get_true_optimal_dose", "simulate_patient_outcome", "run_single_trial_simulation"))

# [Pseudocode] Loop s from 1 to 6 in parallel. Run simulation functions on independent worker environments.
raw_parallel_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
  run_single_trial_simulation(s, compiled_stan_model, N_patients, T_max, target_dlt, 
                              q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi)
})

# Close down the background workers and release memory allocation.
parallel::stopCluster(cl)


# ==============================================================================
# 6. POST-PROCESSING SCORECARD GENERATION
# ==============================================================================
# [Pseudocode] Extract sub-vectors across all list entries inside raw_parallel_results using vector slicing.
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

p_acc_dfcrm  <- sapply(raw_parallel_results, function(x) x$personal_acc_dfcrm)
p_acc_trialr <- sapply(raw_parallel_results, function(x) x$personal_acc_trialr)
p_acc_m5base <- sapply(raw_parallel_results, function(x) x$personal_acc_m5base)
p_acc_m5cov  <- sapply(raw_parallel_results, function(x) x$personal_acc_m5cov)

all_assigned_doses <- unlist(lapply(raw_parallel_results, function(x) x$cov_assigned_doses))
all_patient_risks   <- unlist(lapply(raw_parallel_results, function(x) x$cov_patient_risks))
risk_alignment_r    <- cor(all_patient_risks, all_assigned_doses) # Pearson correlation matrix calculation

# [Human] Bind all demographic logs across the trials to examine the high-risk, low-weight subpopulation safety.
all_demographics <- do.call(rbind, lapply(raw_parallel_results, function(x) x$demographics))
fragile_cohort   <- all_demographics[all_demographics$female == 1 & all_demographics$bmi < 22, ]

cat("\n========================================================================\n")
cat("          UPGRADED CLINICAL METRICS (SMOKE TEST SCORECARD)              \n")
cat("========================================================================\n")

# [Pseudocode] Create table frame, compute across-trial mean scores, round figures to 2 decimals, and print output.
oc_summary <- data.frame(
  Framework             = c("dfcrm (MLE)", "trialr (MCMC)", "M5Base (Pure Pop)", "M5Cov (Personalized)"),
  Global_MTD_Accuracy   = c(mean(final_mtd_dfcrm == 3), mean(final_mtd_trialr == 3), mean(final_mtd_m5base == 3), mean(final_mtd_m5cov == 3)) * 100,
  Personalization_Acc   = c(mean(p_acc_dfcrm), mean(p_acc_trialr), mean(p_acc_m5base), mean(p_acc_m5cov)),
  Avg_DLTs_Per_Trial    = c(mean(total_dlts_dfcrm), mean(total_dlts_trialr), mean(total_dlts_m5base), mean(total_dlts_m5cov)),
  True_Patient_Overdose = c(mean(pt_overdose_dfcrm), mean(pt_overdose_trialr), mean(pt_overdose_m5base), mean(pt_overdose_m5cov))
)
oc_summary[, 2:5] <- round(oc_summary[, 2:5], 2)
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
  # [Human] Re-enable default screen/sleep profiles so your computer doesn't stay awake infinitely after you walk away.
  cat("\nRestoring default Windows sleep timers...\n")
  system("powercfg /change monitor-timeout-ac 15", ignore.stdout = TRUE)
  system("powercfg /change standby-timeout-ac 15", ignore.stdout = TRUE)
}

