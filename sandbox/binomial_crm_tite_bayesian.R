library(rethinking)
library(MASS)
set.seed(123)

# ==============================================================================
# 1. PUBLISHED REALISTIC TRIAL ARCHITECTURE & DATA SETUP
# ==============================================================================
# A 25-patient trial with staggered arrivals and 28-day DLT observation windows
N      <- 25
T_max  <- 28

# Reconstructed patient assignments and DLT outcomes across 5 sequential cohorts
d_assigned <- c(1,1,1,  2,2,2,  3,3,3,  4,4,4,4,  3,3,3,  4,4,4,  5,5,5,5,5,5)
y_observed <- c(0,0,0,  0,0,0,  0,0,1,  1,1,0,0,  0,0,0,  0,1,0,  1,1,0,1,0,0)

# Realized follow-up status at the time of an interim analysis window:
# If y == 1, DLT timestamp is caught early. If y == 0, some are fully cleared (28), 
# while recently accrued patients are partially censored (e.g., 5 to 22 days)
t_followup <- c(28,28,28, 28,28,28, 28,28,14, 8,12,28,28, 28,28,21, 28,19,28, 6,11,28,4,28,15)

# Derive weights: DLT events provide immediate full evidence (1.0)
w_weights  <- ifelse(y_observed == 1, 1.0, t_followup / T_max)

# Generate a continuous patient-level biological exposure metric (e.g., Centered PK AUC)
# Higher dose levels generally lead to higher systemic drug exposure
pk_auc_raw <- rnorm(N, mean = d_assigned * 50, sd = 25)
pk_auc_scaled <- as.numeric(scale(pk_auc_raw))

# Model Priors and target
q       <- c(0.05, 0.12, 0.25, 0.40, 0.55) # Prior Clinical Skeleton
logit_q <- logit(q)                        # Logit anchors
target  <- 0.25

# Pack global data list
trial_data <- list(
  y = y_observed, d = d_assigned, w = w_weights, 
  pk_auc = pk_auc_scaled, prior_means = logit_q
)

# Initialize output tables
comparison_matrix <- matrix(NA, nrow=5, ncol=5)
rownames(comparison_matrix) <- c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5")
colnames(comparison_matrix) <- c("M1_Bin_Logistic", "M2_Classic_CRM", "M3_BWS_Snapshot", "M4_TiTE_BWS", "M5_TiTE_Reg_Cov")

# ==============================================================================
# MODEL 1: STANDARD BINOMIAL LOGISTIC REGRESSION (Snapshot View)
# ==============================================================================
m1_logistic <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a + b * d,
    a ~ dnorm(0, 1.5),
    b ~ dnorm(0, 1)
  ), data = list(y = y_observed, d = d_assigned)
)
samples_m1 <- extract.samples(m1_logistic, n=10000)
p_post_m1  <- sapply(1:5, function(dose) plogis(samples_m1$a + samples_m1$b * dose))
comparison_matrix[, 1] <- colMeans(p_post_m1)

# ==============================================================================
# MODEL 2: CLASSIC ONE-PARAMETER POWER CRM (Snapshot View)
# ==============================================================================
m2_crm <- quap(
  alist(
    y ~ dbinom(1, p),
    p <- exp(log_q * exp(theta)),
    theta ~ dnorm(0, 1.34)
  ), data = list(y = y_observed, log_q = log(q[d_assigned])), start = list(theta = 0)
)
samples_m2 <- extract.samples(m2_crm, n=10000)
p_post_m2  <- sapply(1:5, function(dose) q[dose]^exp(samples_m2$theta))
comparison_matrix[, 2] <- colMeans(p_post_m2)

# ==============================================================================
# MODEL 3: BETA-WEIGHT SKELETON (BWS) GROUP PRIORS (Snapshot View)
# ==============================================================================
m3_bws <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a[d],
    a[1] ~ dnorm(-2.94, 0.2), # logit(0.05)
    a[2] ~ dnorm(-1.99, 0.2), # logit(0.12)
    a[3] ~ dnorm(-1.10, 0.2), # logit(0.25)
    a[4] ~ dnorm(-0.41, 0.2), # logit(0.40)
    a[5] ~ dnorm( 0.20, 0.2)  # logit(0.55)
  ), data = list(y = y_observed, d = d_assigned), start = list(a = logit_q)
)
p_post_m3  <- plogis(extract.samples(m3_bws, depth=2, n=10000)$a)
comparison_matrix[, 3] <- colMeans(p_post_m3)

# ==============================================================================
# MODEL 4: TIME-TO-EVENT BWS (TiTE-BWS via Bounded Optim)
# ==============================================================================
m4_tite_bws_posterior <- function(parameters, data) {
  alpha <- parameters[1:5]
  
  log_p      <- plogis(alpha[data$d], log.p = TRUE)
  log_one_p  <- plogis(alpha[data$d], lower.tail = FALSE, log.p = TRUE)
  
  # Likelihood accounting for right-censored exposure weights
  log_lik    <- sum( data$y * log_p + (1 - data$y) * data$w * log_one_p )
  log_prior  <- sum(dnorm(alpha, mean = data$prior_means, sd = 0.2, log = TRUE))
  
  return(-(log_lik + log_prior))
}

fit_m4 <- optim(par = logit_q, fn = m4_tite_bws_posterior, method = "L-BFGS-B",
                lower = rep(-5, 5), upper = rep(5, 5), hessian = TRUE, data = trial_data)
samples_m4 <- mvrnorm(n=10000, mu=fit_m4$par, Sigma=solve(fit_m4$hessian))
p_post_m4  <- plogis(samples_m4)
comparison_matrix[, 4] <- colMeans(p_post_m4)

# ==============================================================================
# MODEL 5: TiTE-BWS REGRESSION COVARIATE ADJUSTMENT (No Intercept Form)
# ==============================================================================
m5_tite_reg_posterior <- function(parameters, data) {
  alpha   <- parameters[1:5] # Absolute baseline positions per dose
  beta_pk <- parameters[6]   # Continuous biomarker coefficient
  
  # Compute individual patient risk linear combination
  eta_i     <- alpha[data$d] + (beta_pk * data$pk_auc)
  log_p     <- plogis(eta_i, log.p = TRUE)
  log_one_p <- plogis(eta_i, lower.tail = FALSE, log.p = TRUE)
  
  log_lik   <- sum( data$y * log_p + (1 - data$y) * data$w * log_one_p )
  log_prior <- sum(dnorm(alpha, mean = data$prior_means, sd = 0.2, log = TRUE)) + 
    dnorm(beta_pk, mean = 0, sd = 1.0, log = TRUE)
  
  return(-(log_lik + log_prior))
}

fit_m5 <- optim(par = c(logit_q, 0), fn = m5_tite_reg_posterior, method = "L-BFGS-B",
                lower = c(rep(-5, 5), -3), upper = c(rep(5, 5), 3), hessian = TRUE, data = trial_data)
samples_m5 <- mvrnorm(n=10000, mu=fit_m5$par, Sigma=solve(fit_m5$hessian))

# For direct spatial comparison against the others, we isolate the baseline dose effects
p_post_m5_baseline <- plogis(samples_m5[, 1:5])
comparison_matrix[, 5] <- colMeans(p_post_m5_baseline)

# ==============================================================================
# 2. GENERATE COMPARATIVE VISUALIZATION
# ==============================================================================
par(mar=c(5, 4, 4, 10) + 0.1)
plot(1:5, q, type="b", pch=19, lty=3, ylim=c(0, 0.70), lwd=2, col="gray50",
     xlab="Dose Level", ylab="Probability of DLT", main="Evolution of Phase I Trial Models", las=1)

# Overlay individual architectural iterations
lines(1:5, comparison_matrix[,1], col="deepskyblue",   lwd=2, type="o", pch=1)
lines(1:5, comparison_matrix[,2], col="orange",        lwd=2, type="o", pch=2)
lines(1:5, comparison_matrix[,3], col="purple",        lwd=2, type="o", pch=3)
lines(1:5, comparison_matrix[,4], col="forestgreen",   lwd=3, type="o", pch=17)
lines(1:5, comparison_matrix[,5], col="firebrick",     lwd=3, type="o", pch=15)

abline(h=target, lty=2, col="grey30") # Target Toxicity Boundary

legend(x = 5.2, y = 0.7, xpd = TRUE, bty = "n", cex = 0.8,
       legend=c("Prior Skeleton", "1. Binomial Logistic", "2. Power CRM", 
                "3. BWS Snapshot", "4. TiTE-BWS", "5. TiTE-BWS + PK"),
       col=c("gray50", "deepskyblue", "orange", "purple", "forestgreen", "firebrick"),
       lty=c(3,1,1,1,1,1), lwd=c(2,2,2,2,3,3), pch=c(19,1,2,3,17,15))

# ==============================================================================
# 3. INTERIM DECISION READOUT
# ==============================================================================
cat("\n================= POSTERIOR PROBABILITY MATRIX =================")
print(round(comparison_matrix, 3))
cat("================================================================\n\n")

cat("--- MAXIMUM TOLERATED DOSE (MTD) RECOMMENDATIONS ---\n")
cat("1. Binomial Logistic GLM MTD : Dose", which.min(abs(comparison_matrix[,1] - target)), "\n")
cat("2. Classic Power-CRM MTD     : Dose", which.min(abs(comparison_matrix[,2] - target)), "\n")
cat("3. Baseline BWS Snapshot MTD : Dose", which.min(abs(comparison_matrix[,3] - target)), "\n")
cat("4. Corrected TiTE-BWS MTD    : Dose", which.min(abs(comparison_matrix[,4] - target)), "\n")
cat("5. TiTE Regression (PK) MTD  : Dose", which.min(abs(comparison_matrix[,5] - target)), "\n\n")

# Report isolated covariate findings from Model 5
post_beta_pk <- samples_m5[, 6]
cat("--- COVARIATE BIOMARKER INTERPRETATION (Model 5) ---\n")
cat("Beta PK Effect Size (Mean) : ", round(mean(post_beta_pk), 3), "\n")
cat("Beta PK 95% Credible Int.  : [", round(quantile(post_beta_pk, 0.025), 3), ",", round(quantile(post_beta_pk, 0.975), 3), "]\n")

#######################################################
###### adding a stan section for tite regression

library(rstan)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ==============================================================================
# 1. GENERATE GENERALIZED DESIGN MATRIX STRUCTURE
# ==============================================================================
# Coerce the assigned doses into an explicit unordered factor
df_trial <- data.frame(
  y_resp = y_observed,
  fDose  = as.factor(d_assigned), 
  pk_auc = pk_auc_scaled
)

# Crucial Formula: '0 + fDose' strips out the global intercept. 
# This automatically expands the factor into 5 mutually exclusive binary indicator columns.
mModMatrix <- model.matrix(y_resp ~ 0 + fDose + pk_auc, data = df_trial)

# Let's inspect what the design matrix looks like:
# Col 1-5: Indicator flags for Dose 1, 2, 3, 4, 5
# Col 6  : Continuous scaled PK covariate values
print(head(mModMatrix))

# ==============================================================================
# 2. RUN STAN SAMPLING ENGINE
# ==============================================================================
lStanData <- list(
  Ntotal      = length(df_trial$y_resp),
  Ncol        = ncol(mModMatrix),
  X           = mModMatrix,
  y           = df_trial$y_resp,
  w           = w_weights,
  prior_means = logit_q
)

# Compile the external stan file
stanDso <- rstan::stan_model(file = "sandbox/TiTE_likelihood_regression.stan")

# Execute Hamiltonian Monte Carlo sampling
fit.stan <- sampling(
  stanDso, 
  data    = lStanData, 
  iter    = 2000, 
  chains  = 4, 
  cores   = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 13)
)

# ==============================================================================
# 3. INTERPRET POSTERIOR EXTRUCTIONS
# ==============================================================================
# Print out summary statistics directly from the chain tracking
print(fit.stan, pars = c("betas"))

# Extract standard matrix arrays
mCoef <- rstan::extract(fit.stan)$betas
colnames(mCoef) <- colnames(mModMatrix)

# Separate absolute base dose-odds from clinical regression covariates
mDoseBaselines <- mCoef[, 1:5]
mCovariates    <- mCoef[, 6:ncol(mCoef), drop = FALSE]

# Convert the absolute dose baseline posteriors back to standard probability metrics
mPostProbabilities <- plogis(mDoseBaselines)

cat("\n--- STAN HMC MEAN BASELINE DOSES (PROBABILITY SCALE) ---\n")
print(round(colMeans(mPostProbabilities), 3))

cat("\n--- STAN HMC COVARIATE EFFECT INFERENCE ---\n")
for(col in colnames(mCovariates)) {
  samples <- mCovariates[, col]
  cat(col, ": Mean =", round(mean(samples), 3), 
      " | 95% CI = [", round(quantile(samples, 0.025), 3), ",", round(quantile(samples, 0.975), 3), "]\n")
}
