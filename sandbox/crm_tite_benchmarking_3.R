library(dfcrm)
library(trialr)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ==============================================================================
# 1. COMPILE THE MODEL ONCE OUTSIDE THE SIMULATION
# ==============================================================================
cat("Compiling your custom Stan model... (Reused natively)\n")
compiled_stan_model <- rstan::stan_model(file = "sandbox/TiTE_likelihood_regression.stan")

# ==============================================================================
# 2. EXPERIMENTAL & SIMULATION SETUP
# ==============================================================================
set.seed(123)
N_patients <- 20
T_max      <- 28
target_dlt <- 0.25
q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton)

# Hidden True Biological Profile
true_alpha <- qlogis(c(0.04, 0.10, 0.24, 0.42, 0.60)) 
true_b_sex <- 0.8   # Females have higher toxicity
true_b_bmi <- -0.3  # Lower BMI increases toxicity

# Generate Incoming Patient Stream
patients <- data.frame(
  id          = 1:N_patients,
  female      = sample(c(0, 1), N_patients, replace = TRUE),
  bmi_raw     = rnorm(N_patients, mean = 26, sd = 4),
  arrival_day = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE)))
)
patients$bmi_centered <- patients$bmi_raw - 26

# ==============================================================================
# HELPER FUNCTION: DYNAMIC TRUE OUTCOME SIMULATOR
# ==============================================================================
simulate_patient_outcome <- function(dose, female, bmi_centered, arrival_day) {
  logit_p <- true_alpha[dose] + (true_b_sex * female) + (true_b_bmi * bmi_centered)
  p_true  <- plogis(logit_p)
  
  dlt <- rbinom(1, 1, p_true)
  if (dlt == 1) {
    dlt_day <- sample(1:24, 1)
    clearance_day <- arrival_day + dlt_day
  } else {
    clearance_day <- arrival_day + T_max
  }
  return(list(dlt = dlt, clearance_day = clearance_day))
}

# ==============================================================================
# 3. RUN FOUR INDEPENDENT PARALLEL TRIAL REALITIES
# ==============================================================================

# Vectors to track independent histories
d_dfcrm  <- rep(NA, N_patients); y_dfcrm  <- rep(NA, N_patients); c_dfcrm  <- rep(NA, N_patients)
d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients)
d_m5base <- rep(NA, N_patients); y_m5base <- rep(NA, N_patients); c_m5base <- rep(NA, N_patients) # CONTROL
d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients) # FULL COVARIATE

cat("\nExecuting 4-way parallel trial matrix loops...\n")

for (i in 1:N_patients) {
  current_time <- patients$arrival_day[i]
  
  # ----------------------------------------------------------------------------
  # REALITY 1: THE DFCRM WORLD (SOTA Frequentist TiTE-CRM)
  # ----------------------------------------------------------------------------
  if (i == 1) {
    d_dfcrm[i] <- 1
  } else {
    comp <- 1:(i-1)
    days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
    days_on_study <- ifelse(y_dfcrm[comp] == 1, c_dfcrm[comp] - patients$arrival_day[comp], days_on_study)
    w_dfcrm       <- ifelse(y_dfcrm[comp] == 1, 1.0, days_on_study / T_max)
    
    fit_dfcrm <- dfcrm::titecrm(prior=q_skeleton, target=target_dlt, tox=y_dfcrm[comp], 
                                level=d_dfcrm[comp], weights=w_dfcrm, followup=T_max, model="empiric")
    d_dfcrm[i] <- min(fit_dfcrm$mtd, d_dfcrm[i-1] + 1)
  }
  out_dfcrm <- simulate_patient_outcome(d_dfcrm[i], patients$female[i], patients$bmi_centered[i], current_time)
  y_dfcrm[i] <- out_dfcrm$dlt; c_dfcrm[i] <- out_dfcrm$clearance_day
  
  # --------------------------------------------------------------------------
  # REALITY 2: THE TRIALR WORLD (SOTA Bayesian TiTE-CRM via Stan)
  # --------------------------------------------------------------------------
  if (i == 1) {
    d_trialr[i] <- 1
  } else {
    comp <- 1:(i-1)
    days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
    days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
    w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
    
    fit_trialr <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", 
                                   doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=1.34, refresh=0)
    rec_dose <- which.min(abs(fit_trialr$prob_tox - target_dlt))
    d_trialr[i] <- min(rec_dose, d_trialr[i-1] + 1)
  }
  out_trialr <- simulate_patient_outcome(d_trialr[i], patients$female[i], patients$bmi_centered[i], current_time)
  y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day
  
  # ----------------------------------------------------------------------------
  # REALITY 3: M5 BASE (Your Engine, NO Covariates - THE SCIENTIFIC CONTROL GROUP)
  # ----------------------------------------------------------------------------
  if (i == 1) {
    d_m5base[i] <- 1
  } else {
    comp <- 1:(i-1)
    days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
    days_on_study <- ifelse(y_m5base[comp] == 1, c_m5base[comp] - patients$arrival_day[comp], days_on_study)
    w_m5base      <- ifelse(y_m5base[comp] == 1, 1.0, days_on_study / T_max)
    
    df_base <- data.frame(fDose=factor(d_m5base[comp], levels=1:5))
    X_base  <- model.matrix(~ 0 + fDose, data=df_base) # Doses only! Ncol = 5.
    
    stan_data_base <- list(Ntotal=nrow(X_base), Ncol=ncol(X_base), X=X_base, y=as.array(y_m5base[comp]), w=as.array(w_m5base), prior_means=logit_q)
    fit_stan_base  <- rstan::sampling(compiled_stan_model, data=stan_data_base, iter=1500, chains=2, refresh=0, warmup=500)
    
    base_alphas <- rstan::extract(fit_stan_base)$betas
    
    pred_probs_base <- sapply(1:5, function(d) {
      mean(plogis(base_alphas[, d])) # Pure population curves
    })
    
    rec_dose <- which.min(abs(pred_probs_base - target_dlt))
    d_m5base[i] <- min(rec_dose, d_m5base[i-1] + 1)
  }
  out_m5base <- simulate_patient_outcome(d_m5base[i], patients$female[i], patients$bmi_centered[i], current_time)
  y_m5base[i] <- out_m5base$dlt; c_m5base[i] <- out_m5base$clearance_day
  
  # ----------------------------------------------------------------------------
  # REALITY 4: M5 COVARIATE (Your Engine, WITH Sex & BMI Personalization)
  # ----------------------------------------------------------------------------
  if (i == 1) {
    d_m5cov[i] <- 1
  } else {
    comp <- 1:(i-1)
    days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
    days_on_study <- ifelse(y_m5cov[comp] == 1, c_m5cov[comp] - patients$arrival_day[comp], days_on_study)
    w_m5cov       <- ifelse(y_m5cov[comp] == 1, 1.0, days_on_study / T_max)
    
    df_cov <- data.frame(fDose=factor(d_m5cov[comp], levels=1:5), female=patients$female[comp], bmi=patients$bmi_centered[comp])
    X_cov  <- model.matrix(~ 0 + fDose + female + bmi, data=df_cov) # Ncol = 7
    
    stan_data_cov <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov[comp]), w=as.array(w_m5cov), prior_means=logit_q)
    fit_stan_cov  <- rstan::sampling(compiled_stan_model, data=stan_data_cov, iter=1500, chains=2, refresh=0, warmup=500)
    
    cov_betas  <- rstan::extract(fit_stan_cov)$betas
    inc_female <- patients$female[i]; inc_bmi <- patients$bmi_centered[i]
    
    pred_probs_cov <- sapply(1:5, function(d) {
      eta_pred <- cov_betas[, d] + (cov_betas[, 6] * inc_female) + (cov_betas[, 7] * inc_bmi)
      mean(plogis(eta_pred))
    })
    
    rec_dose <- which.min(abs(pred_probs_cov - target_dlt))
    d_m5cov[i]  <- min(rec_dose, d_m5cov[i-1] + 1)
  }
  out_m5cov <- simulate_patient_outcome(d_m5cov[i], patients$female[i], patients$bmi_centered[i], current_time)
  y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day
}

# ==============================================================================
# 4. REPORT TRUE INTERACTIVE COMPARISON MATRIX
# ==============================================================================
results <- data.frame(
  Pt_ID    = patients$id,
  Female   = patients$female,
  BMI      = round(patients$bmi_raw, 1),
  D_dfcrm  = d_dfcrm,  Y_dfcrm  = y_dfcrm,
  D_trialr = d_trialr, Y_trialr = y_trialr,
  D_M5Base = d_m5base, Y_M5Base = y_m5base,  # Control Group
  D_M5Cov  = d_m5cov,  Y_M5Cov  = y_m5cov   # Personalized Group
)
cat("\n================================================================================\n")
cat("          VERIFIED METHODOLOGICAL 4-WAY PARALLEL TRIAL SIMULATION             \n")
cat("================================================================================\n")
print(results)
