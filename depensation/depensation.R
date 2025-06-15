# A depensatory stock recruitment model
# h controls the strength and onset of depensation:
# - It is the spawner biomass where survival drops to half its maximum
# - Small h → depensation occurs only at very low abundance
# - Large h → depensation affects a broader range of spawner levels

# The modified recruitment equations are:
#
#   Ricker:         R(S) = [S / (h + S)] * S * exp(a - b * S)
#   Beverton-Holt:  R(S) = ([S / (h + S)] * S * exp(a)) / (1 + b * S)
#
# This approach allows flexible control over the onset of depensation:
# - Small h → depensation kicks in only at very low spawner abundance
# - Large h → depensation persists over a broader range of S
#
# ln(R/S) is plotted to visualize recruits-per-spawner and the replacement line

ln_a <- 2.5
b <- 0.01
h <- c(10, 50, 200)
col <- c("black", "dodgerblue3", "darkorange")
s <- seq(0.01, 500, length.out = 300) # avoid zero for log scale

## helper functions
ricker <- function(S, h) (S / (h + S)) * S * exp(ln_a - b * S)
bevholt <- function(S, h) ((S / (h + S)) * S * exp(ln_a)) / (1 + b * S)

## matrices: rows = S, cols = different h values
R_ricker <- sapply(h, ricker, S = s)
R_bh <- sapply(h, bevholt, S = s)

lnRPS_ricker <- log(R_ricker / s)
lnRPS_bh <- log(R_bh / s)

## 4-panel plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1) Ricker: R vs S
plot(s, R_ricker[, 1],
  type = "n", ylim = c(0, max(R_ricker)),
  xlab = "Spawning stock biomass (S)", ylab = "Recruitment (R)",
  main = "Ricker: R vs S"
)
for (i in seq_along(h)) lines(s, R_ricker[, i], col = col[i], lwd = 2)
legend("topright", legend = paste("h =", h), col = col, lwd = 2, bty = "n")

# 2) Ricker: ln(R/S) vs S
plot(s, lnRPS_ricker[, 1],
  type = "n", ylim = range(lnRPS_ricker),
  xlab = "S", ylab = "ln(R/S)", main = "Ricker: recruits-per-spawner"
)
abline(h = 0, col = "gray70", lty = 2)
for (i in seq_along(h)) lines(s, lnRPS_ricker[, i], col = col[i], lwd = 2)

# 3) Beverton-Holt: R vs S
plot(s, R_bh[, 1],
  type = "n", ylim = c(0, max(R_bh)),
  xlab = "S", ylab = "Recruitment (R)",
  main = "Beverton-Holt: R vs S"
)
for (i in seq_along(h)) lines(s, R_bh[, i], col = col[i], lwd = 2)
legend("bottomright", legend = paste("h =", h), col = col, lwd = 2, bty = "n")

# 4) Beverton-Holt: ln(R/S) vs S
plot(s, lnRPS_bh[, 1],
  type = "n", ylim = range(lnRPS_bh),
  xlab = "S", ylab = "ln(R/S)", main = "Beverton-Holt: recruits-per-spawner"
)
abline(h = 0, col = "gray70", lty = 2)
for (i in seq_along(h)) lines(s, lnRPS_bh[, i], col = col[i], lwd = 2)

#-------------------------------------------------------------------------------
# Do a simulation and make sure it works
set.seed(123)

## true values -------------------------------------------------------------
ln_a_true <- 0.25 # ln-alpha
b_true <- 0.01 # density dependence
h_true <- 50 # half-saturation for survival term
sigma <- 0.35 # sd of log error

n_obs <- 100 # observations per data set
n_sim <- 1000 # monte-carlo replicates

## storage for mles --------------------------------------------------------
ln_a_hat <- numeric(n_sim)
b_hat <- numeric(n_sim)
h_hat <- numeric(n_sim)

## negative log-likelihood -------------------------------------------------
nll <- function(par, S, R) {
  ln_a <- par[1]
  b <- par[2]
  h <- exp(par[3]) # positivity with exp
  sd <- exp(par[4])
  mu <- log(S) + log(S / (h + S)) + ln_a - b * S
  -sum(dnorm(log(R), mu, sd, log = TRUE))
}

## monte-carlo loop --------------------------------------------------------
for (k in 1:n_sim) {
  S <- runif(n_obs, 1, 500)
  mu <- log(S) + log(S / (h_true + S)) + ln_a_true - b_true * S
  R <- exp(rnorm(n_obs, mu, sigma))

  start <- c(ln_a = 0, b = 0.005, ln_h = log(30), ln_sd = log(0.3))
  fit <- optim(start, nll, S = S, R = R, method = "BFGS")
  ln_a_hat[k] <- fit$par[1]
  b_hat[k] <- fit$par[2]
  h_hat[k] <- exp(fit$par[3])

  if (k %% 100 == 0) cat("finished", k, "of", n_sim, "\n")
}

## plots -------------------------------------------------------------------
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

hist(ln_a_hat,
  breaks = 20, col = "gray85",
  main = expression("mle of a"), xlab = expression(hat(a))
)
abline(v = ln_a_true, col = "dodgerblue3", lwd = 3, lty = 2)

hist(b_hat,
  breaks = 20, col = "gray85",
  main = expression("mle of b"), xlab = expression(hat(b))
)
abline(v = b_true, col = "dodgerblue3", lwd = 3, lty = 2)

hist(h_hat,
  breaks = 20, col = "gray85",
  main = expression("mle of" ~ h), xlab = expression(hat(h))
)
abline(v = h_true, col = "dodgerblue3", lwd = 3, lty = 2)
