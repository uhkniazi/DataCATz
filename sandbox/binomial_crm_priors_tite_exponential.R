library(rethinking)
library(MASS) # Required for sampling from multivariate normal distributions
set.seed(42)

# ==============================================================================
# 1. SIMULATE TRIAL DATA WITH TIME TRACKING
# ==============================================================================
N <- 50
d <- rep(1:5, each = 10)                  # 50 patients total across 5 doses
p_true <- c(0.05, 0.10, 0.20, 0.30, 0.45) # True, unknown full-period probability
q <- c(0.05, 0.10, 0.20, 0.35, 0.50)       # Clinical Skeleton Prior
logit_q <- logit(q)                        # Logit-skeleton anchors

# Define the maximum observation window (e.g., 28 days)
T_max <- 28

# Draw a true, latent day-of-toxicity using an exponential distribution based on p_true
# True continuous hazard rate lambda derived from p_true: p = 1 - exp(-lambda * T_max)
lambda_true <- -log(1 - p_true) / T_max
t_toxicity <- rexp(N, rate = lambda_true[d])

# Simulate random enrollment times / monitoring observation cutoffs
# Some patients were just enrolled (short tracking), others are finished (full tracking)
t_observed_available <- runif(N, min = 1, max = T_max)

# Realized outcome: DLT occurs only if true toxicity day happens BEFORE their observed tracking window
y <- ifelse(t_toxicity <= t_observed_available, 1, 0)

# Exact time recorded for the trial analysis:
# If DLT happened, we stop the clock on that day. If not, it's their current censored follow-up time.
t_trial <- ifelse(y == 1, t_toxicity, t_observed_available)

# Calculate the fractional exposure weights for the TiTE model
# If DLT occurred: weight is 1.0 immediately. If not: it's the fraction of time navigated.
w <- ifelse(y == 1, 1.0, t_trial / T_max)

# Group standard datasets
dat <- list(y = y, d = as.numeric(d))
dat_crm <- list(y = y, d = d, log_q = log(q[d]))
dat_bws <- list(y = y, d = d)

# ==========================================
# MODEL 1: Standard Logistic Regression (Aggregated Data)
# ==========================================
m_logistic <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a + b * d,
    a ~ dcauchy(0, 2.5),
    b ~ dcauchy(0, 2.5)
  ), data = dat
)
p_post_logistic <- sapply(1:5, function(dose) plogis(extract.samples(m_logistic)$a + extract.samples(m_logistic)$b * dose))

# ==========================================
# MODEL 2: Classic CRM Power Model (Aggregated Data)
# ==========================================
m_crm <- quap(
  alist(
    y ~ dbinom(1, p),
    p <- exp(log_q * exp(theta)),
    theta ~ dnorm(0, 1)
  ), data = dat_crm, start = list(theta = 0)
)
p_post_crm <- sapply(1:5, function(dose) q[dose]^exp(extract.samples(m_crm)$theta))

# ==========================================
# MODEL 3: Original BWS Group Prior Model (Aggregated Data)
# ==========================================
m_bws <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a[d],
    a[1] ~ dnorm(-2.94, 0.2), 
    a[2] ~ dnorm(-2.20, 0.2), 
    a[3] ~ dnorm(-1.38, 0.2), 
    a[4] ~ dnorm(-0.62, 0.2), 
    a[5] ~ dnorm(0.00, 0.2)   
  ), data = dat_bws, start = list(a = logit_q)
)
p_post_bws <- plogis(extract.samples(m_bws, depth=2)$a)

# ==========================================
# MODEL 4: CUSTOM TiTE-BWS (Weighted Likelihood via Optim)
# ==========================================
tite_bws_posterior <- function(parameters, data) {
  a <- parameters 
  p_i <- 1 / (1 + exp(-a[data$d]))
  
  # Corrected TiTE Likelihood: Time scales evidence via log1p(-p) stability
  log_lik <- sum( data$y * log(p_i + 1e-15) + (1 - data$y) * data$w * log1p(-p_i) )
  log_prior <- sum(dnorm(a, mean = data$prior_means, sd = 0.2, log = TRUE))
  return(-(log_lik + log_prior)) # Negative because optim minimizes
}

fit_tite <- optim(par = logit_q, fn = tite_bws_posterior, hessian = TRUE, 
                  data = list(y=y, d=d, w=w, prior_means=logit_q))

# Sample from Laplace posterior covariance matrix
samples_tite <- mvrnorm(n=10000, mu=fit_tite$par, Sigma=solve(fit_tite$hessian))
p_post_tite <- 1 / (1 + exp(-samples_tite))

# ==========================================
# MODEL 5: CONTINUOUS EXPONENTIAL SURVIVAL MODEL
# ==========================================
survival_posterior <- function(parameters, data) {
  a <- parameters
  lambda_i <- exp(a[data$d])
  
  # Continuous Exponential Survival log-likelihood: y*log(hazard) - hazard*time
  log_lik <- sum( data$y * log(lambda_i + 1e-15) - (lambda_i * data$t) )
  log_prior <- sum(dnorm(a, mean = data$prior_means, sd = 0.2, log = TRUE))
  return(-(log_lik + log_prior))
}

# Shift the prior skeleton to a continuous scale anchor: lambda = -log(1-q)/T_max
prior_lambdas_logit <- log(-log(1 - q) / T_max)

fit_surv <- optim(par = prior_lambdas_logit, fn = survival_posterior, hessian = TRUE,
                  data = list(y=y, d=d, t=t_trial, prior_means=prior_lambdas_logit))

samples_surv <- mvrnorm(n=10000, mu=fit_surv$par, Sigma=solve(fit_surv$hessian))
# Project continuous hazards back up to a cumulative 4-week window probability: p = 1 - exp(-lambda * T_max)
p_post_surv <- 1 - exp(-exp(samples_surv) * T_max)

# ==============================================================================
# 5. VISUALIZATION & COMPARATIVE FIGURE
# ==============================================================================
# Set wider margins for layout
par(mar=c(5, 4, 4, 8) + 0.1)

plot(1:5, p_true, type="b", pch=19, lty=2, ylim=c(0, 0.60), cex=1.5,
     xlab="Dose Level", ylab="Toxicity Probability (Over 4 Weeks)", 
     main="Dose-Response Evolution: Snapshot vs. Survival Models", las=1)

# Plot snapshots
lines(1:5, colMeans(p_post_logistic), col="deepskyblue", lwd=2, type="o", pch=1)
lines(1:5, colMeans(p_post_crm), col="firebrick1", lwd=2, type="o", pch=2)
lines(1:5, colMeans(p_post_bws), col="purple", lwd=2, type="o", pch=3)

# Plot Time-to-Event and true survival architectures
lines(1:5, colMeans(p_post_tite), col="darkgreen", lwd=3, type="o", pch=17)
lines(1:5, colMeans(p_post_surv), col="orange", lwd=3, type="o", pch=15)

# Coordinate layout for side legend allocation
legend(x = 5.2, y = 0.6, xpd = TRUE,
       legend=c("True Target Curve", "1. Logistic snapshot", "2. CRM snapshot", 
                "3. BWS snapshot", "4. TiTE-BWS (Weighted)", "5. Exponential Survival"),
       col=c("black", "deepskyblue", "firebrick1", "purple", "darkgreen", "orange"), 
       lwd=c(2,2,2,2,3,3), lty=c(2,1,1,1,1,1), pch=c(19,1,2,3,17,15), bty="n", cex=0.85)

# ==========================================
# FINAL DECISION-MAKING PARSER
# ==========================================
target <- 0.25
cat("\n--- MAXIMUM TOLERATED DOSE (MTD) RECOMMENDATIONS ---\n")
cat("Standard Logistic GLM MTD : Dose", which.min(abs(colMeans(p_post_logistic) - target)), "\n")
cat("Classic Power-CRM MTD     : Dose", which.min(abs(colMeans(p_post_crm) - target)), "\n")
cat("Your Baseline BWS MTD     : Dose", which.min(abs(colMeans(p_post_bws) - target)), "\n")
cat("Corrected TiTE-BWS MTD    : Dose", which.min(abs(colMeans(p_post_tite) - target)), "\n")
cat("Continuous Survival MTD   : Dose", which.min(abs(colMeans(p_post_surv) - target)), "\n")
