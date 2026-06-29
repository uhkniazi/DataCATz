# Name: generate_section_4_toxicity_curves.R
# Auth: u.niazi@soton.ac.uk
# Date: 23/06/2026
# Desc: Generates Section 4 Figure using Base R. Evaluates average fitted curves 
#       under Scenario C across 20 trials for trialr vs M5Cov.
#       Implements Andrew Gelman style contrast coding (-0.5 / +0.5) for sex.

library(dfcrm)
library(trialr)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)

# Compile custom Stan engine once on master
cat("Compiling custom Monotonic Stan engine...\n")
compiled_stan_model <- rstan::stan_model(file = "TiTE_regression_monotonic_skeleton.stan")

# ==============================================================================
# GLOBAL EXPERIMENTAL CONFIGURATION
# ==============================================================================
set.seed(42)
N_sims     <- 20   # Compute the average fitted curve across 20 datasets
N_patients <- 20   
T_max      <- 28   
target_dlt <- 0.25 

q_skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q    <- qlogis(q_skeleton) 

# Scenario C Truth
true_alpha <- qlogis(c(0.02, 0.15, 0.35, 0.55, 0.70)) 
true_b_sex <- 0.8  # Log-odds ratio between female (+0.5) and male (-0.5)
true_b_bmi <- -0.3 

get_true_risk <- function(dose, sex_centered, bmi_centered, true_alpha, true_b_sex, true_b_bmi) {
  plogis(true_alpha[dose] + (true_b_sex * sex_centered) + (true_b_bmi * bmi_centered))
}

get_true_optimal_dose <- function(sex_centered, bmi_centered, true_alpha, true_b_sex, true_b_bmi, target_dlt = 0.25) {
  true_risks <- sapply(1:5, function(d) get_true_risk(d, sex_centered, bmi_centered, true_alpha, true_b_sex, true_b_bmi))
  return(which.min(abs(true_risks - target_dlt)))
}

simulate_patient_outcome <- function(dose, sex_centered, bmi_centered, arrival_day, T_max, true_alpha, true_b_sex, true_b_bmi) {
  p_true  <- get_true_risk(dose, sex_centered, bmi_centered, true_alpha, true_b_sex, true_b_bmi)
  dlt     <- rbinom(1, 1, p_true) 
  clearance_day <- arrival_day + ifelse(dlt == 1, sample(1:24, 1), T_max)
  return(list(dlt = dlt, clearance_day = clearance_day, p_true = p_true))
}

# ==============================================================================
# SIMULATION LOOP (Tracks fitted curves at Patient 20 using Contrast Centering)
# ==============================================================================
run_single_curve_simulation <- function(sim_id, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi) {
  
  set.seed(42 + sim_id)
  patients <- data.frame(
    id           = 1:N_patients,
    # GELMAN-STYLE CONTRAST CODING: -0.5 = Male, +0.5 = Female
    sex_centered = sample(c(-0.5, 0.5), N_patients, replace = TRUE),
    bmi_raw      = rnorm(N_patients, mean = 26, sd = 4),
    arrival_day  = cumsum(c(0, sample(4:10, N_patients - 1, replace = TRUE)))
  )
  patients$bmi_centered <- patients$bmi_raw - 26
  patients$true_optimal <- sapply(1:N_patients, function(idx) {
    get_true_optimal_dose(patients$sex_centered[idx], patients$bmi_centered[idx], true_alpha, true_b_sex, true_b_bmi, target_dlt)
  })
  
  d_trialr <- rep(NA, N_patients); y_trialr <- rep(NA, N_patients); c_trialr <- rep(NA, N_patients)
  d_m5cov  <- rep(NA, N_patients); y_m5cov  <- rep(NA, N_patients); c_m5cov  <- rep(NA, N_patients)
  
  for (i in 1:N_patients) {
    current_time <- patients$arrival_day[i]
    comp <- 1:(i-1)
    
    # 1. trialr Allocation Loop
    if (i == 1) { d_trialr[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_trialr[comp] == 1, c_trialr[comp] - patients$arrival_day[comp], days_on_study)
      w_trialr      <- ifelse(y_trialr[comp] == 1, 1.0, days_on_study / T_max)
      fit_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr[comp], tox=y_trialr[comp], weights=w_trialr, beta_sd=sqrt(1.34), chains=4, cores=1, refresh=0)
      d_trialr[i]   <- min(which.min(abs(fit_trialr$prob_tox - target_dlt)), d_trialr[i-1] + 1)
    }
    out_trialr <- simulate_patient_outcome(d_trialr[i], patients$sex_centered[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_trialr[i] <- out_trialr$dlt; c_trialr[i] <- out_trialr$clearance_day
    
    # 2. M5Cov Allocation Loop
    if (i == 1) { d_m5cov[i] <- 1 } else {
      days_on_study <- pmin(current_time - patients$arrival_day[comp], T_max)
      days_on_study <- ifelse(y_m5cov[comp] == 1, c_m5cov[comp] - patients$arrival_day[comp], days_on_study)
      w_m5cov        <- ifelse(y_m5cov[comp] == 1, 1.0, days_on_study / T_max)
      X_cov         <- model.matrix(~ 0 + factor(d_m5cov[comp], levels=1:5) + sex_centered + bmi_centered, data=patients[comp,])
      stan_data_c   <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov[comp]), w=as.array(w_m5cov), prior_means=logit_q)
      fit_stan_c    <- rstan::sampling(compiled_stan_model, data=stan_data_c, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
      extracted_c   <- rstan::extract(fit_stan_c)
      cov_alpha     <- extracted_c$alpha     
      cov_beta_cov  <- extracted_c$beta_cov  
      pred_c        <- sapply(1:5, function(d) {
        mean(plogis(cov_alpha[, d] + (cov_beta_cov[, 1] * patients$sex_centered[i]) + (cov_beta_cov[, 2] * patients$bmi_centered[i])))
      })
      d_m5cov[i]    <- min(which.min(abs(pred_c - target_dlt)), d_m5cov[i-1] + 1)
    }
    out_m5cov <- simulate_patient_outcome(d_m5cov[i], patients$sex_centered[i], patients$bmi_centered[i], current_time, T_max, true_alpha, true_b_sex, true_b_bmi)
    y_m5cov[i] <- out_m5cov$dlt; c_m5cov[i] <- out_m5cov$clearance_day
  }
  
  # --- FINAL EXTRACTION AT END OF TRIAL (Patient 20) ---
  final_time_context <- patients$arrival_day[N_patients]
  
  # Re-fit trialr one last time to capture final dataset curve
  days_on_study_t <- pmin(final_time_context + T_max - patients$arrival_day, T_max)
  days_on_study_t <- ifelse(y_trialr == 1, c_trialr - patients$arrival_day, days_on_study_t)
  w_trialr_final  <- ifelse(y_trialr == 1, 1.0, days_on_study_t / T_max)
  final_trialr    <- trialr::stan_crm(skeleton=q_skeleton, target=target_dlt, model="empiric", doses_given=d_trialr, tox=y_trialr, weights=w_trialr_final, beta_sd=sqrt(1.34), chains=4, cores=1, refresh=0)
  curve_trialr    <- final_trialr$prob_tox
  
  # Explicitly calculate full-length weights for M5Cov final fit
  days_on_study_c <- pmin(final_time_context + T_max - patients$arrival_day, T_max)
  days_on_study_c <- ifelse(y_m5cov == 1, c_m5cov - patients$arrival_day, days_on_study_c)
  w_m5cov_final   <- ifelse(y_m5cov == 1, 1.0, days_on_study_c / T_max)
  
  X_cov <- model.matrix(~ 0 + factor(d_m5cov, levels=1:5) + sex_centered + bmi_centered, data=patients)
  stan_data_c <- list(Ntotal=nrow(X_cov), Ncol=ncol(X_cov), X=X_cov, y=as.array(y_m5cov), w=as.array(w_m5cov_final), prior_means=logit_q)
  final_m5cov <- rstan::sampling(compiled_stan_model, data=stan_data_c, iter=1000, chains=4, cores=1, refresh=0, warmup=400)
  ext <- rstan::extract(final_m5cov)
  
  # Extract posterior fitted curves under new coordinate definitions
  # Evaluating at (0, 0) yields the true Unisex Midpoint patient profile
  curve_m5cov_base      <- sapply(1:5, function(d) mean(plogis(ext$alpha[, d] + ext$beta_cov[,1]*0 + ext$beta_cov[,2]*0)))
  curve_m5cov_fragile   <- sapply(1:5, function(d) mean(plogis(ext$alpha[, d] + ext$beta_cov[,1]*0.5 + ext$beta_cov[,2]*(-5))))
  curve_m5cov_resilient <- sapply(1:5, function(d) mean(plogis(ext$alpha[, d] + ext$beta_cov[,1]*(-0.5) + ext$beta_cov[,2]*5)))
  
  return(list(
    trialr = curve_trialr,
    m5_base = curve_m5cov_base,
    m5_fragile = curve_m5cov_fragile,
    m5_resilient = curve_m5cov_resilient
  ))
}

# ==============================================================================
# RUN PARALLEL EXECUTIONS & COLLECT CURVES
# ==============================================================================
num_workers <- min(detectCores() - 4, N_sims)
cl <- parallel::makeCluster(num_workers)
parallel::clusterEvalQ(cl, { library(dfcrm); library(trialr); library(rstan) })
parallel::clusterExport(cl, c("compiled_stan_model", "N_patients", "T_max", "target_dlt", 
                              "q_skeleton", "logit_q", "true_alpha", "true_b_sex", "true_b_bmi", 
                              "get_true_risk", "get_true_optimal_dose", "simulate_patient_outcome", "run_single_curve_simulation"))

curve_results <- parallel::parLapply(cl, 1:N_sims, function(s) {
  run_single_curve_simulation(s, compiled_stan_model, N_patients, T_max, target_dlt, q_skeleton, logit_q, true_alpha, true_b_sex, true_b_bmi)
})
parallel::stopCluster(cl)

# Aggregation (Compute means across the 20 simulated iterations with NA handling)
mean_trialr       <- rowMeans(sapply(curve_results, function(x) x$trialr), na.rm = TRUE)
mean_m5_base      <- rowMeans(sapply(curve_results, function(x) x$m5_base), na.rm = TRUE)
mean_m5_fragile   <- rowMeans(sapply(curve_results, function(x) x$m5_fragile), na.rm = TRUE)
mean_m5_resilient <- rowMeans(sapply(curve_results, function(x) x$m5_resilient), na.rm = TRUE)

# Calculate True Theoretical Curves for Comparison Profiles under New Center Points
true_pop       <- sapply(1:5, function(d) get_true_risk(d, 0, 0, true_alpha, true_b_sex, true_b_bmi))
true_fragile   <- sapply(1:5, function(d) get_true_risk(d, 0.5, -5, true_alpha, true_b_sex, true_b_bmi))   # Female (+0.5)
true_resilient <- sapply(1:5, function(d) get_true_risk(d, -0.5, 5, true_alpha, true_b_sex, true_b_bmi))   # Male (-0.5)

# ==============================================================================
# HIGH-RESOLUTION MULTI-PAGE PDF VISUALIZATION ENGINE (FINALIZED)
# ==============================================================================
# Open the PDF graphics device with optimal dimensions for a single-panel layout
pdf("curves.pdf", width = 8, height = 6.5)

# Set up a single-panel layout per page with clean margins
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# --- PAGE 1 / PANEL A: POPULATION CURVE BENDING ---
plot(1:5, q_skeleton, type = "b", pch = 1, lty = 3, col = "gray40", lwd = 2,
     ylim = c(0, 1), xlim = c(1, 5), xaxt = "n", xlab = "Dose Level", 
     ylab = "Probability of Toxicity (DLT)", 
     main = "A: Curve Adaptation Under Misspecification\n(Scenario C: Misspecified Skeleton)", 
     cex.main = 1.2, cex.lab = 1.1)
axis(1, at = 1:5, labels = paste("Dose", 1:5))
abline(h = target_dlt, col = "darkgreen", lty = 2, lwd = 1.2) 

lines(1:5, true_pop, type = "b", pch = 16, col = "black", lwd = 2.5)
lines(1:5, mean_trialr, type = "b", pch = 17, col = "firebrick3", lwd = 2, lty = 4)
lines(1:5, mean_m5_base, type = "b", pch = 15, col = "dodgerblue4", lwd = 2, lty = 1)

# FIXED: Label cleanly updated to show it represents a true balanced Unisex Midpoint
legend("topleft", 
       legend = c("Initial Skeleton Input", "True Population Reality", 
                  "trialr (TiTE-CRM) Mean Fit", "M5Cov Baseline Fit (Midpoint Unisex Profile)"),
       col = c("gray40", "black", "firebrick3", "dodgerblue4"), 
       pch = c(1, 16, 17, 15), lty = c(3, 1, 4, 1), lwd = 2, bty = "n", cex = 0.95)


# --- PAGE 2 / PANEL B: PHENOTYPE PROFILE DIVERGENCE ---
plot(1:5, true_fragile, type = "l", col = "deeppink3", lwd = 2, lty = 2,
     ylim = c(0, 1), xlim = c(1, 5), xaxt = "n", xlab = "Dose Level", 
     ylab = "Probability of Toxicity (DLT)", 
     main = "B: Phenotypic Profile Resolution\n(Scenario C: Misspecified Skeleton)", 
     cex.main = 1.2, cex.lab = 1.1)
axis(1, at = 1:5, labels = paste("Dose", 1:5))
abline(h = target_dlt, col = "darkgreen", lty = 2, lwd = 1.2)

lines(1:5, true_resilient, type = "l", col = "forestgreen", lwd = 2, lty = 2)
lines(1:5, mean_m5_fragile, type = "b", pch = 17, col = "deeppink4", lwd = 2.5)
lines(1:5, mean_m5_resilient, type = "b", pch = 15, col = "darkgreen", lwd = 2.5)
lines(1:5, mean_trialr, type = "l", col = "firebrick3", lwd = 2, lty = 4)

legend("topleft", 
       legend = c("True Fragile (Low-BMI Female)", "M5Cov Fragile Estimate", 
                  "True Resilient (High-BMI Male)", "M5Cov Resilient Estimate",
                  "trialr Blind Average Fit"),
       col = c("deeppink3", "deeppink4", "forestgreen", "darkgreen", "firebrick3"), 
       pch = c(NA, 17, NA, 15, NA), lty = c(2, 1, 2, 1, 4), lwd = 2, bty = "n", cex = 0.95)

# Close and write the PDF file cleanly
dev.off()
cat("Successfully generated contrast-centered two-page vector graphics file: 'curves.pdf'\n")
