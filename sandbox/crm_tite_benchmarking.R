# ==============================================================================
# DEPENDENCY CHECK & INITIALIZATION
# ==============================================================================
# Ensure required production packages are installed:
# install.packages(c("rethinking", "MASS", "rstan", "dfcrm", "trialr"))

library(rethinking)
library(MASS)
library(rstan)
library(dfcrm)
library(trialr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

# ==============================================================================
# 1. THE PHASE I BENCHMARK DATASET SETUP
# ==============================================================================
N     <- 25
T_max <- 28

# Cohort tracking vectors
d_assigned <- c(1,1,1,  2,2,2,  3,3,3,  4,4,4,4,  3,3,3,  4,4,4,  5,5,5,5,5,5)
y_observed <- c(0,0,0,  0,0,0,  0,0,1,  1,1,0,0,  0,0,0,  0,1,0,  1,1,0,1,0,0)
t_followup <- c(28,28,28, 28,28,28, 28,28,14, 8,12,28,28, 28,28,21, 28,19,28, 6,11,28,4,28,15)

# Derive fractional follow-up weights
w_weights  <- ifelse(y_observed == 1, 1.0, t_followup / T_max)

# Standardize patient-level PK covariate
pk_auc_raw    <- rnorm(N, mean = d_assigned * 50, sd = 25)
pk_auc_scaled <- as.numeric(scale(pk_auc_raw))

# Target and Prior Skeleton Profile
q       <- c(0.05, 0.12, 0.25, 0.40, 0.55)
logit_q <- logit(q)
target  <- 0.25

trial_data <- list(
  y = y_observed, d = d_assigned, w = w_weights, 
  pk_auc = pk_auc_scaled, prior_means = logit_q
)

# Structure our master benchmark comparison grid
benchmark_matrix <- matrix(NA, nrow = 5, ncol = 8)
rownames(benchmark_matrix) <- c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5")
colnames(benchmark_matrix) <- c(
  "M2_Custom_CRM", "SOTA_dfcrm_Classic", 
  "M4_Custom_TiTE", "SOTA_dfcrm_TiTE", 
  "M5_Custom_Optim", "M5_Custom_Stan", "SOTA_trialr_Stan",
  "Skeleton"
)
benchmark_matrix[, "Skeleton"] <- q

# ==============================================================================
# 2. RUNNING CUSTOM OPTIM MODELS (M2 & M4)
# ==============================================================================
# --- Model 2: Custom Power CRM ---
m2_crm <- quap(
  alist(
    y ~ dbinom(1, p),
    p <- exp(log_q * exp(theta)),
    theta ~ dnorm(0, 1.34)
  ), data = list(y = y_observed, log_q = log(q[d_assigned])), start = list(theta = 0)
)
benchmark_matrix[, "M2_Custom_CRM"] <- colMeans(sapply(1:5, function(d) q[d]^exp(extract.samples(m2_crm, n=5000)$theta)))

# --- Model 4: Custom TiTE-BWS ---
m4_tite_bws_posterior <- function(parameters, data) {
  alpha <- parameters[1:5]
  log_p      <- plogis(alpha[data$d], log.p = TRUE)
  log_one_p  <- plogis(alpha[data$d], lower.tail = FALSE, log.p = TRUE)
  log_lik    <- sum(data$y * log_p + (1 - data$y) * data$w * log_one_p)
  log_prior  <- sum(dnorm(alpha, mean = data$prior_means, sd = 0.2, log = TRUE))
  return(-(log_lik + log_prior))
}
fit_m4 <- optim(par = logit_q, fn = m4_tite_bws_posterior, method = "L-BFGS-B",
                lower = rep(-5, 5), upper = rep(5, 5), hessian = TRUE, data = trial_data)
benchmark_matrix[, "M4_Custom_TiTE"] <- plogis(fit_m4$par)

# --- Model 5: Custom TiTE-BWS Regression (optim variation) ---
m5_tite_reg_posterior <- function(parameters, data) {
  alpha   <- parameters[1:5]
  beta_pk <- parameters[6]
  eta_i   <- alpha[data$d] + (beta_pk * data$pk_auc)
  log_p     <- plogis(eta_i, log.p = TRUE)
  log_one_p <- plogis(eta_i, lower.tail = FALSE, log.p = TRUE)
  log_lik   <- sum(data$y * log_p + (1 - data$y) * data$w * log_one_p)
  log_prior <- sum(dnorm(alpha, mean = data$prior_means, sd = 0.2, log = TRUE)) + dnorm(beta_pk, mean = 0, sd = 2.5, log = TRUE)
  return(-(log_lik + log_prior))
}
fit_m5 <- optim(par = c(logit_q, 0), fn = m5_tite_reg_posterior, method = "L-BFGS-B",
                lower = c(rep(-5, 5), -5), upper = c(rep(5, 5), 5), hessian = TRUE, data = trial_data)
benchmark_matrix[, "M5_Custom_Optim"] <- plogis(fit_m5$par[1:5])

# ==============================================================================
# 3. RUNNING CUSTOM INLINE STAN ENGINE (M5 Custom MCMC Variation)
# ==============================================================================
stan_model_code <- "
data {
  int<lower=1> Ntotal;
  int<lower=1> Ncol;
  matrix[Ntotal, Ncol] X;
  array[Ntotal] int<lower=0, upper=1> y;
  vector[Ntotal] w;
  vector[5] prior_means;
}
parameters {
  vector[Ncol] betas;
}
model {
  for (d in 1:5) {
    betas[d] ~ normal(prior_means[d], 0.2);
  }
  if (Ncol > 5) {
    for (c in 6:Ncol) {
      betas[c] ~ normal(0, 2.5);
    }
  }
  vector[Ntotal] eta = X * betas;
  for (i in 1:Ntotal) {
    if (y[i] == 1) {
      target += log_inv_logit(eta[i]);
    } else {
      target += w[i] * log1m_inv_logit(eta[i]);
    }
  }
}
"

# Construct our unconstrained cell-means matrix
df_trial   <- data.frame(y_resp = y_observed, fDose = as.factor(d_assigned), pk_auc = pk_auc_scaled)
mModMatrix <- model.matrix(y_resp ~ 0 + fDose + pk_auc, data = df_trial)

lStanData  <- list(Ntotal = nrow(mModMatrix), Ncol = ncol(mModMatrix), X = mModMatrix, y = df_trial$y_resp, w = w_weights, prior_means = logit_q)

# Compile inline directly
stanDso  <- rstan::stan_model(model_code = stan_model_code)
fit.stan <- rstan::sampling(stanDso, data = lStanData, iter = 2000, chains = 4, warmup = 1000, verbose = FALSE, refresh = 0)

mCoef <- rstan::extract(fit.stan)$betas
benchmark_matrix[, "M5_Custom_Stan"] <- colMeans(plogis(mCoef[, 1:5]))

# ==============================================================================
# 4. EXECUTING STATE-OF-THE-ART (SOTA) BENCHMARKS
# ==============================================================================
# SOTA 1: dfcrm - Snapshot Power Model (FIXED: "empirical")
sota_dfcrm_classic <- dfcrm::crm(
  prior  = q, 
  target = target, 
  tox    = y_observed, 
  level  = d_assigned, 
  model  = "empiric"
)
benchmark_matrix[, "SOTA_dfcrm_Classic"] <- sota_dfcrm_classic$ptox

# SOTA 2: dfcrm - Time-to-Event Model (FIXED: "empirical")
sota_dfcrm_tite <- dfcrm::titecrm(
  prior    = q, 
  target   = target, 
  tox      = y_observed, 
  level    = d_assigned, 
  weights  = w_weights,  # Pass our pre-calculated TiTE weights directly!
  followup = T_max, 
  model    = "empiric"
)
benchmark_matrix[, "SOTA_dfcrm_TiTE"] <- sota_dfcrm_tite$ptox

# SOTA 3: trialr - Stan Bayesian Empirical CRM
sota_trialr <- trialr::stan_crm(
  skeleton    = q, 
  target      = target, 
  model       = "empiric", 
  doses_given = d_assigned, 
  tox         = y_observed, 
  weights     = w_weights, 
  beta_sd     = 1.34,
  refresh     = 0
)

benchmark_matrix[, "SOTA_trialr_Stan"] <- sota_trialr$prob_tox

target_dlt = target
# ==============================================================================
# 5. FINAL SIDE-BY-SIDE BENCHMARK READOUT
# ==============================================================================
cat("\n============================================================================\n")
cat("                  FINAL PHASE I MODEL BENCHMARK MATRIX                       \n")
cat("============================================================================\n")
print(round(benchmark_matrix, 3))
cat("============================================================================\n\n")

cat("--- MAXIMUM TOLERATED DOSE (MTD) DETERMINATION ---\n")
cat("M2 Custom CRM MTD        : Dose", which.min(abs(benchmark_matrix[, "M2_Custom_CRM"] - target_dlt)), "\n")
cat("SOTA dfcrm Classic MTD   : Dose", sota_dfcrm_classic$mtd, "\n")
cat("M4 Custom TiTE MTD       : Dose", which.min(abs(benchmark_matrix[, "M4_Custom_TiTE"] - target_dlt)), "\n")
cat("SOTA dfcrm TiTE MTD      : Dose", sota_dfcrm_tite$mtd, "\n")
cat("M5 Custom Stan MTD       : Dose", which.min(abs(benchmark_matrix[, "M5_Custom_Stan"] - target_dlt)), "\n")
cat("SOTA trialr MTD          : Dose", which.min(abs(benchmark_matrix[, "SOTA_trialr_Stan"] - target_dlt)), "\n\n")

# Report Covariate extraction matching
cat("--- ISOLATED COVARIATE CHECK (MODEL 5 VS STANDARDS) ---\n")
cat("Custom Stan Beta PK (Mean) : ", round(mean(mCoef[, 6]), 3), "\n")
cat("Custom Stan PK 95% CI      : [", round(quantile(mCoef[, 6], 0.025), 3), ",", round(quantile(mCoef[, 6], 0.975), 3), "]\n")
cat("Note: Industry packages like 'dfcrm' or base 'trialr::stan_crm' do not support\n")
cat("direct patient-level covariate design-matrix embedding natively without script hacking.\n")