# binomial_and_crm.R

library(rethinking)
set.seed(42)

# dose levels (categorical)
d <- rep(1:5, each = 10)   # 50 patients total

# "true" toxicity probabilities (unknown in real life)
p_true <- c(0.05, 0.10, 0.20, 0.30, 0.45)

# generate outcomes
y <- rbinom(length(d), size = 1, prob = p_true[d])

table(d, y)

dat <- list(
  y = y,
  d = as.numeric(d)
)

m_logistic <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a + b * d,
    a ~ dcauchy(0, 2.5),
    b ~ dcauchy(0, 2.5)
  ),
  data = dat
)

d_seq <- 1:5
post <- extract.samples(m_logistic, n = 1000)

p_post_logistic <- sapply(d_seq, function(dose) {
  plogis(post$a + post$b * dose)
})

apply(p_post_logistic, 2, mean)

### crm
## These are prior guesses of toxicity probabilities
q <- c(0.05, 0.10, 0.20, 0.35, 0.50)

dat_crm <- list(
  y = y,
  d = d,
  q = q
)

# m_crm <- quap(
#   alist(
#     y ~ dbinom(1, p),
#     p <- q[d]^(exp(theta)),
#     theta ~ dnorm(0, 1)
#   ),
#   data = dat_crm,
#   start = list(theta = 0)
# )

# # 1. Define the custom log-posterior function
# crm_log_posterior <- function(theta) {
#   # Retrieve data from environment
#   y <- dat_crm$y
#   d <- dat_crm$d
#   q <- dat_crm$q
#   
#   # The model equations
#   p <- q[d]^(exp(theta))
#   
#   # Log-Likelihood + Log-Prior
#   log_lik <- sum(dbinom(y, size = 1, prob = p, log = TRUE))
#   log_prior <- dnorm(theta, mean = 0, sd = 1, log = TRUE)
#   
#   return(log_lik + log_prior)
# }
# 
# # 2. Pass the function directly into quap
# m_crm <- quap(crm_log_posterior, start = list(theta = 0))
# 
# # 3. View your summary
# precis(m_crm)

# 1. Pre-calculate the log of the skeleton values for each patient's dose
# This removes the problematic q[d] index lookup from the formula
dat_crm$log_q <- log(dat_crm$q[dat_crm$d])

# 2. Rewrite the formula using pure algebra
# p = q^exp(theta)  ->  log(-log(p)) = theta + log(-log(q))
# Or simply calculate p using the pre-computed log_q variable:
m_crm <- quap(
  alist(
    y ~ dbinom(1, p),
    p <- exp(log_q * exp(theta)),
    theta ~ dnorm(0, 1)
  ),
  data = dat_crm,
  start = list(theta = 0)
)

precis(m_crm)

post_crm <- extract.samples(m_crm)

p_post_crm <- sapply(1:5, function(dose) {
  q[dose]^exp(post_crm$theta)
})

apply(p_post_crm, 2, mean)


## compare 2 models
plot(1:5, p_true, type="b", pch=16, ylim=c(0,0.6),
     xlab="Dose", ylab="Toxicity prob")

lines(1:5, colMeans(p_post_logistic), col="blue", lwd=2)
lines(1:5, colMeans(p_post_crm), col="red", lwd=2)

legend("topleft",
       legend=c("True", "Logistic", "CRM"),
       col=c("black", "blue", "red"),
       lwd=2)


### decision step (what trial actually does)
target <- 0.25

abs_diff <- abs(colMeans(p_post_crm) - target)
which.min(abs_diff)


#############################################################
###### second version of models with a hierarchical binomial model
library(rethinking)
set.seed(42)

# 1. Simulate Trial Data
d <- rep(1:5, each = 10)   # 50 patients total
p_true <- c(0.05, 0.10, 0.20, 0.30, 0.45) # Unknown in real life
## toxicities at different doses, based on the probabilities
y <- rbinom(length(d), size = 1, prob = p_true[d])

# 2. Setup Data Lists
dat <- list(y = y, d = as.numeric(d))
q <- c(0.05, 0.10, 0.20, 0.35, 0.50) # Clinical Skeleton Prior

# ==========================================
# MODEL 1: Standard Logistic Regression
# ==========================================
m_logistic <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a + b * d,
    a ~ dcauchy(0, 2.5),
    b ~ dcauchy(0, 2.5)
  ), data = dat
)
## posterior predictive simulation
p_post_logistic <- sapply(1:5, function(dose) plogis(extract.samples(m_logistic)$a + extract.samples(m_logistic)$b * dose))

# ==========================================
# MODEL 2: Classic CRM (Power Model)
# ==========================================
dat_crm <- list(y = y, d = d, log_q = log(q[d]))
m_crm <- quap(
  alist(
    y ~ dbinom(1, p),
    p <- exp(log_q * exp(theta)),
    theta ~ dnorm(0, 1)
  ), data = dat_crm, start = list(theta = 0)
)
post_crm <- extract.samples(m_crm)
p_post_crm <- sapply(1:5, function(dose) q[dose]^exp(post_crm$theta))

# ==========================================
# MODEL 3: YOUR IDEA (Beta-Weight Skeleton as Group Priors)
# ==========================================
# Convert skeleton probabilities to logit scale coordinates
logit_q <- logit(q) 

dat_bws <- list(y = y, d = d)

m_bws <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a[d],
    # 5 independent group priors tightly locked around the skeleton positions
    a[1] ~ dnorm(-2.94, 0.2), # logit(0.05)
    a[2] ~ dnorm(-2.20, 0.2), # logit(0.10)
    a[3] ~ dnorm(-1.38, 0.2), # logit(0.20)
    a[4] ~ dnorm(-0.62, 0.2), # logit(0.35)
    a[5] ~ dnorm(0.00, 0.2)   # logit(0.50)
  ), data = dat_bws, start = list(a = logit_q)
)
p_post_bws <- plogis(extract.samples(m_bws, depth=2)$a)

# ==========================================
# VISUALIZATION & COMPARISON
# ==========================================
plot(1:5, p_true, type="b", pch=16, ylim=c(0,0.6), xlab="Dose", ylab="Toxicity prob", main="Model Comparison vs Truth")
lines(1:5, colMeans(p_post_logistic), col="blue", lwd=2)
lines(1:5, colMeans(p_post_crm), col="red", lwd=2)
lines(1:5, colMeans(p_post_bws), col="darkgreen", lwd=2) # Your model

legend("topleft", legend=c("True", "Logistic", "CRM", "Your Model (BWS)"),
       col=c("black", "blue", "red", "darkgreen"), lwd=2)

# ==========================================
# FINAL TRIAL DECISIONS
# ==========================================
target <- 0.25
cat("CRM MTD Recommendation:", which.min(abs(colMeans(p_post_crm) - target)), "\n")
cat("Your Model MTD Recommendation:", which.min(abs(colMeans(p_post_bws) - target)), "\n")
