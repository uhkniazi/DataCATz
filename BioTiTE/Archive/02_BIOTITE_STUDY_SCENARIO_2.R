# ==============================================================================
# BIOTITE STUDY 1: SCENARIO 2 (CONFOUNDED DOSE-BIOMARKER ARCHITECTURE)
# Auth: u.niazi@soton.ac.uk
# Date: 08/09/2026
# Desc: Evaluates Tier 0 (GLM), Tier 1 (Standard Stan Normal(0,2)), and Tier 2 
#       (BioTiTE) when Gene1 expression is partially confounded with dose level.
# ==============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

# ------------------------------------------------------------------------------
# 1. PARAMETER SETUP & SCENARIO 2 DATA GENERATION
# ------------------------------------------------------------------------------
n       <- 40
n_doses <- 5
n_genes <- 5

# Target Skeleton & Logit Intercepts
skeleton   <- c(0.05, 0.10, 0.20, 0.35, 0.50)
alpha_true <- qlogis(skeleton)

# Ground Truth Gene Effects (2 Causal, 3 Null)
beta_true <- c(
  Gene1 = 1.0,  # Strong causal (Confounding target)
  Gene2 = 0.5,  # Moderate causal
  Gene3 = 0.0,  # Null
  Gene4 = 0.0,  # Null
  Gene5 = 0.0   # Null
)

# Phase I Dose Allocation (Skewed to lower/mid doses)
dose_probs <- c(0.30, 0.30, 0.20, 0.10, 0.10)
dose       <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)

# --- SCENARIO 2: Dose-Dependent Gene1 Prevalence ---
# Gene1 probability increases with assigned dose cohort
g1_prob <- c(0.10, 0.20, 0.40, 0.60, 0.80)

G <- matrix(rnorm(n * n_genes), nrow = n, ncol = n_genes)
G[, 1] <- rbinom(n, size = 1, prob = g1_prob[dose]) # Confounded binary trait
colnames(G) <- paste0("Gene", 1:n_genes)

# Linear Predictor & DLT Outcome Generation
eta <- alpha_true[dose] + as.vector(G %*% beta_true)
p   <- plogis(eta)
dlt <- rbinom(n, size = 1, prob = p)

df_sim <- data.frame(
  patient_id = 1:n,
  dose       = dose,
  G,
  p_true     = round(p, 3),
  dlt        = dlt
)

cat("=== SCENARIO 2 DATA SUMMARY ===\n")
cat("Dose Distribution & Observed DLTs:\n")
print(table(Dose = df_sim$dose, DLT = df_sim$dlt))

cat("\nGene 1 Prevalence by Assigned Dose Level:\n")
print(round(tapply(df_sim$Gene1, df_sim$dose, mean), 2))


# ------------------------------------------------------------------------------
# 2. DESIGN MATRIX & STAN DATA ASSEMBLY
# ------------------------------------------------------------------------------
X_design <- model.matrix(~ 0 + factor(dose, levels = 1:5) + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, data = df_sim)

stan_data_biotite <- list(
  Ntotal      = nrow(X_design),
  Ncol        = ncol(X_design),
  X           = X_design,
  y           = as.array(df_sim$dlt),
  w           = rep(1.0, nrow(df_sim)),
  prior_means = alpha_true
)

stan_data_std <- list(
  Ntotal = nrow(X_design),
  Ncol   = ncol(X_design),
  X      = X_design,
  y      = as.array(df_sim$dlt)
)


# ------------------------------------------------------------------------------
# 3. TIER 0: CLASSICAL MAXIMUM LIKELIHOOD REGRESSION (R glm)
# ------------------------------------------------------------------------------
fit_glm <- glm(
  dlt ~ 0 + factor(dose) + Gene1 + Gene2 + Gene3 + Gene4 + Gene5,
  data   = df_sim,
  family = binomial(link = "logit")
)
glm_coef <- summary(fit_glm)$coefficients


# ------------------------------------------------------------------------------
# 4. TIER 1: STANDARD BAYESIAN LOGISTIC MODEL (Normal(0, 2) Regularized)
# ------------------------------------------------------------------------------
compiled_std <- rstan::stan_model(file='binomialRegression.stan')

fit_std <- rstan::sampling(
  compiled_std,
  data    = stan_data_std,
  iter    = 2000,
  warmup  = 1000,
  chains  = 4,
  refresh = 0
)


# ------------------------------------------------------------------------------
# 5. TIER 2: BIOTITE RETROSPECTIVE MONOTONIC MODEL
# ------------------------------------------------------------------------------
compiled_biotite <- rstan::stan_model(file = "TiTE_regression_monotonic_skeleton.stan")

fit_biotite <- rstan::sampling(
  compiled_biotite,
  data    = stan_data_biotite,
  iter    = 2000,
  warmup  = 1000,
  chains  = 4,
  refresh = 0
)


# ------------------------------------------------------------------------------
# 6. COMPARATIVE EXTRACTION & DIAGNOSTIC SUMMARY
# ------------------------------------------------------------------------------
ext_std     <- rstan::extract(fit_std)
ext_biotite <- rstan::extract(fit_biotite)

# --- Gene Coefficient Comparison ---
std_gene_means <- apply(ext_std$betas[, 6:10], 2, mean)
std_gene_sds   <- apply(ext_std$betas[, 6:10], 2, sd)

bio_gene_means <- apply(ext_biotite$beta_cov, 2, mean)
bio_gene_sds   <- apply(ext_biotite$beta_cov, 2, sd)

summary_genes <- data.frame(
  True_Effect   = beta_true,
  GLM_Est       = round(glm_coef[6:10, "Estimate"], 3),
  GLM_SE        = round(glm_coef[6:10, "Std. Error"], 3),
  StdStan_Mean  = round(std_gene_means, 3),
  StdStan_SD    = round(std_gene_sds, 3),
  BioTiTE_Mean  = round(bio_gene_means, 3),
  BioTiTE_SD    = round(bio_gene_sds, 3)
)
rownames(summary_genes) <- paste0("Gene", 1:n_genes)

# --- Dose Intercept & Probability Scale Comparison ---
std_dose_probs <- apply(plogis(ext_std$betas[, 1:5]), 2, mean)
bio_dose_probs <- apply(plogis(ext_biotite$alpha), 2, mean)

summary_dose <- data.frame(
  True_Skeleton = skeleton,
  GLM_Prob      = round(plogis(glm_coef[1:5, "Estimate"]), 3),
  StdStan_Prob  = round(std_dose_probs, 3),
  BioTiTE_Prob  = round(bio_dose_probs, 3)
)
rownames(summary_dose) <- paste0("Dose_", 1:n_doses)

# --- Print Comparative Tables ---
cat("\n======================================================================\n")
cat("SCENARIO 2 RESULTS: DOSE-RESPONSE MONOTONICITY & RISK ESTIMATION\n")
cat("======================================================================\n")
print(summary_dose)

cat("\n======================================================================\n")
cat("SCENARIO 2 RESULTS: BIOMARKER RECOVERY UNDER DOSE-GENE CONFOUNDING\n")
cat("======================================================================\n")
print(summary_genes)
