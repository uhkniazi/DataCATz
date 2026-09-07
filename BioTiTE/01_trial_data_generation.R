# Name: 01_trial_data_generation.R
# Auth: u.niazi@soton.ac.uk
# Date: 07/09/2026
# Desc: TiTE-CRM Models used in Phase 1 clinical trials to assess dosage based
#       on toxicity skeleton. Trial data generation scenario


# ==============================================================================
# BIOTITE STUDY 1: SINGLE DEBUG DATASET GENERATOR
# ==============================================================================
set.seed(42)

n <- 40
n_doses <- 5
n_genes <- 5

# 1. Target Skeleton and True Dose Intercepts
skeleton   <- c(0.05, 0.10, 0.20, 0.35, 0.50)
alpha_true <- qlogis(skeleton)

# 2. True Gene Effects (2 Causal, 3 Null)
beta_true <- c(
  g1 = 1.0,  # Strong causal
  g2 = 0.5,  # Moderate causal
  g3 = 0.0,  # Null
  g4 = 0.0,  # Null
  g5 = 0.0   # Null
)

# 3. Realistic Phase I Ascending Dose Allocation
# Skewed toward lower/mid doses rather than uniform sampling
dose_probs <- c(0.30, 0.30, 0.20, 0.10, 0.10)
dose       <- sample(1:n_doses, size = n, replace = TRUE, prob = dose_probs)

# 4. Standardized Genomic Expression Matrix (N x 5)
# Drawn from standard normal distribution
G <- matrix(rnorm(n * n_genes), nrow = n, ncol = n_genes)
colnames(G) <- paste0("Gene", 1:n_genes)

# 5. Linear Predictor & Toxicity Generation
# eta_i = alpha_{d_i} + sum(beta_j * G_{ij})
eta <- alpha_true[dose] + as.vector(G %*% beta_true)
p   <- plogis(eta)
dlt <- rbinom(n, size = 1, prob = p)

# 6. Assemble Final Cohort Data Frame
df_sim <- data.frame(
  patient_id = 1:n,
  dose       = dose,
  G,
  p_true     = round(p, 3),
  dlt        = dlt
)

# Print Summary Table
cat("=== SIMULATED DOSE DISTRIBUTION & OBSERVED DLTS ===\n")
print(table(Dose = df_sim$dose, DLT = df_sim$dlt))

# # ==============================================================================
# # BASE R DIAGNOSTIC VISUALIZATIONS
# # ==============================================================================
# par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
# 
# # Panel A: Observed DLT Rate vs. Assigned Dose
# obs_rates <- aggregate(dlt ~ dose, data = df_sim, FUN = mean)
# counts    <- table(df_sim$dose)
# 
# plot(obs_rates$dose, obs_rates$dlt, type = "b", pch = 16, col = "firebrick3", lwd = 2,
#      ylim = c(0, 1), xlim = c(1, 5), xaxt = "n",
#      xlab = "Dose Level", ylab = "Observed DLT Rate",
#      main = "A: Observed DLT Rate vs Dose")
# axis(1, at = 1:5, labels = paste0("D1\n(n=", as.vector(counts), ")"))
# abline(h = 0.20, lty = 2, col = "gray50")
# 
# # Panel B: Target Skeleton Curve vs Logit Intercept Reality
# plot(1:5, skeleton, type = "b", pch = 17, col = "dodgerblue4", lwd = 2, lty = 2,
#      ylim = c(0, 1), xaxt = "n",
#      xlab = "Dose Level", ylab = "Probability",
#      main = "B: Target Clinical Skeleton")
# axis(1, at = 1:5, labels = paste("Dose", 1:5))
# points(obs_rates$dose, obs_rates$dlt, pch = 16, col = "firebrick3")
# legend("topleft", legend = c("Target Skeleton", "Simulated Empirical"),
#        col = c("dodgerblue4", "firebrick3"), pch = c(17, 16), lty = c(2, 1), bty = "n", cex = 0.9)
# 
# # Panel C: Observed DLT Outcome vs Primary Causal Feature (Gene 1)
# plot(df_sim$Gene1, df_sim$dlt, pch = 16, col = ifelse(df_sim$dlt == 1, "firebrick3", "black"),
#      xlab = "Gene 1 Expression (Std Dev)", ylab = "DLT Status (0/1)",
#      main = "C: Toxicity vs Gene 1 Signal")
# curve(plogis(alpha_true[3] + beta_true[1] * x), add = TRUE, col = "darkgreen", lwd = 2)
# legend("topleft", legend = c("No DLT (0)", "DLT (1)", "Dose 3 Effect Curve"),
#        col = c("black", "firebrick3", "darkgreen"), pch = c(16, 16, NA), lty = c(NA, NA, 1), bty = "n", cex = 0.85)
# 
# par(mfrow = c(1, 1))

# ==============================================================================
# UPDATED DIAGNOSTIC VISUALIZATIONS (Multi-Gene & Multi-Dose Overlay)
# ==============================================================================
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

# --- Panel A: Observed DLT Rate vs Assigned Dose ---
obs_rates <- aggregate(dlt ~ dose, data = df_sim, FUN = mean)
counts    <- table(factor(df_sim$dose, levels = 1:5)) # Ensures all 5 doses represented

plot(obs_rates$dose, obs_rates$dlt, type = "b", pch = 16, col = "firebrick3", lwd = 2,
     ylim = c(0, 1), xlim = c(1, 5), xaxt = "n",
     xlab = "Dose Level", ylab = "Observed DLT Rate",
     main = "A: DLT Rate vs Assigned Dose")
axis(1, at = 1:5, labels = paste0("D", 1:5, "\n(n=", as.vector(counts), ")"))
abline(h = 0.20, lty = 2, col = "gray50")

# --- Panel B: Target Skeleton vs Simulated Empirical Points ---
plot(1:5, skeleton, type = "b", pch = 17, col = "dodgerblue4", lwd = 2, lty = 2,
     ylim = c(0, 1), xaxt = "n",
     xlab = "Dose Level", ylab = "Probability",
     main = "B: Target Skeleton vs Empirical")
axis(1, at = 1:5, labels = paste("Dose", 1:5))
points(obs_rates$dose, obs_rates$dlt, pch = 16, col = "firebrick3")
legend("topleft", legend = c("Target Skeleton", "Empirical DLTs"),
       col = c("dodgerblue4", "firebrick3"), pch = c(17, 16), lty = c(2, 1), bty = "n", cex = 0.85)

# --- Panel C: Empty Spacer for Layout Symmetry ---
plot.new()

# --- Panel D: Gene 1 (Strong Effect: Beta = 1.0) ---
plot(df_sim$Gene1, df_sim$dlt, pch = 16, col = ifelse(df_sim$dlt == 1, "firebrick3", "black"),
     xlab = "Gene 1 Expression (Std Dev)", ylab = "DLT Status (0/1)",
     main = "D: Gene 1 (Strong: Beta = 1.0)")
# Overlay Dose 1 (Low), Dose 3 (Mid), Dose 5 (High) reference curves
curve(plogis(alpha_true[1] + beta_true[1] * x), add = TRUE, col = "forestgreen", lwd = 2, lty = 3)
curve(plogis(alpha_true[3] + beta_true[1] * x), add = TRUE, col = "darkgreen", lwd = 2)
curve(plogis(alpha_true[5] + beta_true[1] * x), add = TRUE, col = "firebrick4", lwd = 2, lty = 2)
legend("topleft", legend = c("Dose 1 Curve", "Dose 3 Curve", "Dose 5 Curve"),
       col = c("forestgreen", "darkgreen", "firebrick4"), lty = c(3, 1, 2), lwd = 1.5, bty = "n", cex = 0.75)

# --- Panel E: Gene 2 (Weak Effect: Beta = 0.5) ---
plot(df_sim$Gene2, df_sim$dlt, pch = 16, col = ifelse(df_sim$dlt == 1, "firebrick3", "black"),
     xlab = "Gene 2 Expression (Std Dev)", ylab = "DLT Status (0/1)",
     main = "E: Gene 2 (Weak: Beta = 0.5)")
curve(plogis(alpha_true[1] + beta_true[2] * x), add = TRUE, col = "forestgreen", lwd = 2, lty = 3)
curve(plogis(alpha_true[3] + beta_true[2] * x), add = TRUE, col = "darkgreen", lwd = 2)
curve(plogis(alpha_true[5] + beta_true[2] * x), add = TRUE, col = "firebrick4", lwd = 2, lty = 2)

# --- Panel F: Gene 3 (Null Effect: Beta = 0.0) ---
plot(df_sim$Gene3, df_sim$dlt, pch = 16, col = ifelse(df_sim$dlt == 1, "firebrick3", "black"),
     xlab = "Gene 3 Expression (Std Dev)", ylab = "DLT Status (0/1)",
     main = "F: Gene 3 (Null: Beta = 0.0)")
curve(plogis(alpha_true[1] + beta_true[3] * x), add = TRUE, col = "forestgreen", lwd = 2, lty = 3)
curve(plogis(alpha_true[3] + beta_true[3] * x), add = TRUE, col = "darkgreen", lwd = 2)
curve(plogis(alpha_true[5] + beta_true[3] * x), add = TRUE, col = "firebrick4", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# ==============================================================================
# STAN DATA STRUCTURE ASSEMBLY (BioTiTE Retrospective)
# ==============================================================================

# Construct full design matrix: 5 dose indicators + 5 gene expression columns
X_design <- model.matrix(~ 0 + factor(dose, levels = 1:5) + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, data = df_sim)

stan_data <- list(
  Ntotal      = nrow(X_design),
  Ncol        = ncol(X_design),
  X           = X_design,
  y           = as.array(df_sim$dlt),
  w           = rep(1.0, nrow(df_sim)), # Retrospective complete follow-up (w_i = 1)
  prior_means = alpha_true              # Logit-scale target skeleton
)

# Verification Check
cat("\nStan Data Package Ready:\n")
cat("N subjects:", stan_data$Ntotal, "\n")
cat("Design matrix dimensions:", dim(stan_data$X), "\n")


# # --- 3. STAN DATA LIST ASSEMBLY ---
# # Design matrix: 5 dose indicators + 5 gene columns
# X_mat <- model.matrix(~ 0 + factor(dose, levels = 1:5) + Gene1 + Gene2 + Gene3 + Gene4 + Gene5, data = df_sim)
# 
# stan_data <- list(
#   Ntotal      = nrow(X_mat),
#   Ncol        = ncol(X_mat),
#   X           = X_mat,
#   y           = df_sim$dlt,
#   w           = rep(1.0, n), # Static retrospective follow-up
#   prior_means = q_target
# )

# --- 4. MODEL EXECUTION ---
# Compiles and samples from your TiTE_regression_monotonic_skeleton.stan file
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
compiled_model <- rstan::stan_model(file = "TiTE_regression_monotonic_skeleton.stan")

fit <- rstan::sampling(
  compiled_model,
  data    = stan_data,
  iter    = 2000,
  warmup  = 1000,
  chains  = 4,
  refresh = 0
)

print(fit, digits=3)
traceplot(fit, pars='alpha')
traceplot(fit, pars='beta_cov')
# --- 5. POSTERIOR RECOVERY DIAGNOSTICS ---
ext <- rstan::extract(fit)

# Summarize Gene Coefficients (beta_cov)
beta_summary <- t(apply(ext$beta_cov, 2, function(x) {
  c(Mean = mean(x), SD = sd(x), `2.5%` = quantile(x, 0.025), `97.5%` = quantile(x, 0.975))
}))
rownames(beta_summary) <- paste0("Gene", 1:n_genes)

cat("\n=== GENE COEFFICIENT RECOVERY (Target: G1=1.0, G2=0.5, G3-5=0.0) ===\n")
print(round(beta_summary, 3))

# Summarize Dose Intercepts (alpha)
alpha_summary <- t(apply(ext$alpha, 2, function(x) {
  c(Posterior_Logit = mean(x), Fitted_Prob = mean(plogis(x)), True_Skeleton = 0)
}))
rownames(alpha_summary) <- paste0("Dose_", 1:n_doses)
alpha_summary[, "True_Skeleton"] <- skeleton

cat("\n=== DOSE INTERCEPT RECOVERY (Target Skeleton Probabilities) ===\n")
print(round(alpha_summary, 3))
