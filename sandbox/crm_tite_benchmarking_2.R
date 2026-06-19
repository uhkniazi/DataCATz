library(dfcrm)
library(trialr)
library(rstan)
library(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ==============================================================================
# 1. COMPILE YOUR CUSTOM STAN MODEL ONCE (REUSE NATIVELY)
# ==============================================================================
cat("Compiling your custom Stan model... (This happens only once)\n")
compiled_stan_model <- rstan::stan_model(file = "sandbox/TiTE_likelihood_regression.stan")

# ==============================================================================
# 2. DEFINE TRUE OPERATING SCENARIOS & CLINICAL SETUP
# ==============================================================================
set.seed(42)
N_patients <- 24
T_max      <- 28
target_dlt <- 0.25
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton)

# Hidden True Biological Profile (Dose effects if patient is an average Male)
true_alpha <- qlogis(c(0.04, 0.10, 0.24, 0.42, 0.60)) 
true_b_sex <- 0.8   # Females experience higher toxicity
true_b_bmi <- -0.3  # Lower BMI experiences higher toxicity

# ==============================================================================
# 3. GENERATE A COHORT OF VIRTUAL PATIENTS WITH TRAITS
# ==============================================================================
sim_demographics <- data.frame(
  patient_id = 1:N_patients,
  female     = sample(c(0, 1), N_patients, replace = TRUE),
  bmi_raw    = rnorm(N_patients, mean = 26, sd = 4),
  arrival_day= cumsum(c(0, sample(5:15, N_patients - 1, replace = TRUE)))
)
sim_demographics$bmi_centered <- sim_demographics$bmi_raw - 26

# Track dynamic trial trajectories (Driven by your Personalized Model 5 decisions)
assigned_doses <- rep(NA, N_patients)
dlt_outcomes   <- rep(NA, N_patients)
clearance_days <- rep(NA, N_patients)

# Data matrices to log what SOTA packages *would* have chosen at each decision point
dfcrm_recommendations   <- rep(NA, N_patients)
trialr_recommendations  <- rep(NA, N_patients)

# ==============================================================================
# 4. RUN DYNAMIC SEQUENTIAL ENROLLMENT LOOP
# ==============================================================================
cat("\nStarting active clinical simulation loop...\n")

for (i in 1:N_patients) {
  current_time <- sim_demographics$arrival_day[i]
  
  if (i == 1) {
    next_dose <- 1 
    dfcrm_recommendations[i]  <- 1
    trialr_recommendations[i] <- 1
  } else {
    completed_patients <- 1:(i-1)
    y_snapshot <- dlt_outcomes[completed_patients]
    d_snapshot <- assigned_doses[completed_patients]
    
    # Calculate TiTE operational tracking times
    days_on_study <- current_time - sim_demographics$arrival_day[completed_patients]
    days_on_study <- pmin(days_on_study, T_max)
    days_on_study <- ifelse(y_snapshot == 1, clearance_days[completed_patients] - sim_demographics$arrival_day[completed_patients], days_on_study)
    w_snapshot    <- ifelse(y_snapshot == 1, 1.0, days_on_study / T_max)
    
    # --------------------------------------------------------------------------
    # COMPARISON ENGINE A: SOTA dfcrm (Population-level MLE TiTE)
    # --------------------------------------------------------------------------
    fit_dfcrm <- dfcrm::titecrm(
      prior    = q_skeleton, target = target_dlt, 
      tox      = y_snapshot, level = d_snapshot, 
      weights  = w_snapshot, followup = T_max, model = "empiric"
    )
    dfcrm_recommendations[i] <- fit_dfcrm$mtd
    
    # --------------------------------------------------------------------------
    # COMPARISON ENGINE B: SOTA trialr (Population-level Stan MCMC CRM)
    # --------------------------------------------------------------------------
    # trialr::stan_crm expects standard CRM formatting. Since base trialr doesn't natively 
    # handle dynamic TiTE partial weights per row without hacking, we pass the current outcomes.
    fit_trialr <- trialr::stan_crm(
      skeleton    = q_skeleton, 
      target      = target_dlt, 
      model       = "empiric", 
      doses_given = d_snapshot, 
      tox         = y_snapshot, 
      beta_sd     = 1.34, 
      refresh     = 0
    )
    # Extract trialr MTD recommendation based on its population posterior mean
    trialr_recommendations[i] <- which.min(abs(fit_trialr$prob_tox - target_dlt))
    
    # --------------------------------------------------------------------------
    # ENGINE C: Your Custom Personalized Model 5 (Covariates + TiTE via Stan)
    # --------------------------------------------------------------------------
    df_snapshot <- data.frame(
      fDose  = factor(d_snapshot, levels = 1:5),
      female = sim_demographics$female[completed_patients],
      bmi    = sim_demographics$bmi_centered[completed_patients]
    )
    X_snapshot <- model.matrix(~ 0 + fDose + female + bmi, data = df_snapshot)
    
    stan_data <- list(
      Ntotal      = nrow(X_snapshot),
      Ncol        = ncol(X_snapshot),
      X           = X_snapshot,
      y           = as.array(y_snapshot),  # FORCE 1-D ARRAY
      w           = as.array(w_snapshot),  # FORCE 1-D ARRAY
      prior_means = logit_q
    )
    
    fit_stan <- rstan::sampling(
      compiled_stan_model, data = stan_data, 
      iter = 1500, chains = 2, refresh = 0, warmup = 500
    )
    
    stan_betas      <- rstan::extract(fit_stan)$betas
    incoming_female <- sim_demographics$female[i]
    incoming_bmi    <- sim_demographics$bmi_centered[i]
    
    # Calculate personalized toxicity probability for patient i's exact body profile
    pred_probs_incoming <- sapply(1:5, function(d) {
      eta_pred <- stan_betas[, d] + (stan_betas[, 6] * incoming_female) + (stan_betas[, 7] * incoming_bmi)
      mean(plogis(eta_pred))
    })
    
    next_dose <- which.min(abs(pred_probs_incoming - target_dlt))
  }
  
  # ----------------------------------------------------------------------------
  # EXECUTE PERSONALIZED PATIENT ASSIGNMENT & SIMULATE OUTCOME
  # ----------------------------------------------------------------------------
  assigned_doses[i] <- next_dose
  
  logit_p_true <- true_alpha[next_dose] + 
    (true_b_sex * sim_demographics$female[i]) + 
    (true_b_bmi * sim_demographics$bmi_centered[i])
  p_true            <- plogis(logit_p_true)
  dlt_outcomes[i]   <- rbinom(1, size = 1, prob = p_true)
  clearance_days[i] <- sim_demographics$arrival_day[i] + ifelse(dlt_outcomes[i] == 1, sample(1:24, 1), T_max)
}

# ==============================================================================
# 5. GENERATE FINAL BENCHMARK DATA REPORT
# ==============================================================================
sim_demographics$assigned_dose_M5 <- assigned_doses
sim_demographics$rec_dfcrm_TiTE   <- dfcrm_recommendations
sim_demographics$rec_trialr_Stan   <- trialr_recommendations
sim_demographics$dlt_observed     <- dlt_outcomes

cat("\n================================================================\n")
cat("            COMPLETE THREE-WAY FRAMEWORK BENCHMARK              \n")
cat("================================================================\n")
print(round(sim_demographics[, c("patient_id", "female", "bmi_raw", "rec_dfcrm_TiTE", "rec_trialr_Stan", "assigned_dose_M5", "dlt_observed")], 2))
