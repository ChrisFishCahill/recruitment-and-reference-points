# life history inputs
vul <- c(0.05, 0.174, 0.571, 0.899, 1.000, 0.6)
wa <- c(0.05, 0.104, 0.145, 0.186, 0.212, 0.253)
mat <- c(0.05, 0.20, 0.82, 0.99, 1.00, 1.00)
M <- 0.4
n_ages <- length(mat)

# survivorship
lo <- numeric(n_ages)
lo[1] <- 1
for (i in 2:n_ages) lo[i] <- lo[i - 1] * exp(-M)
lo[n_ages] <- lo[n_ages - 1] / (1 - exp(-M))

# recruitment params
Ro <- 10
recK <- 5
sbro <- sum(lo * mat * wa)
ln_alpha <- log(recK / sbro)
br <- (ln_alpha + log(sbro)) / (Ro * sbro)
ninit <- Ro * lo

# simulation settings
n_years <- 100
sdr <- 0.6
set.seed(199)
wt <- rnorm(n_years - 1, 0, sdr)

# objective function using Ft = cslope * (vb - lrp)/vb
get_yield_hcr <- function(par, wt, upow) {
  lrp <- par[1]
  cslope <- par[2]

  log_n <- matrix(-Inf, nrow = n_years, ncol = n_ages)
  log_n[1, ] <- log(ninit)
  ssb <- yield <- vul_bio <- numeric(n_years)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    Ft <- max(0, cslope * (vul_bio[t - 1] - lrp) / vul_bio[t - 1])
    Zt <- Ft * vul + M
    log_n[t, 1] <- ln_alpha + log(ssb[t - 1]) - br * ssb[t - 1] + wt[t - 1]
    for (a in 2:n_ages) {
      log_n[t, a] <- log_n[t - 1, a - 1] - Zt[a - 1]
    }
    log_n[t, n_ages] <- log(
      exp(log_n[t, n_ages]) +
        exp(log_n[t - 1, n_ages]) * exp(-Zt[n_ages])
    )
    ssb[t] <- sum(exp(log_n[t, ]) * mat * wa)
    vul_bio[t] <- sum(exp(log_n[t, ]) * wa * vul)
    yield[t] <- sum(exp(log_n[t, ]) * wa * Ft * vul / Zt * (1 - exp(-Zt)))
  }
  -mean((yield[-1])^upow)
}

# thinking abou hyperbolic absolute risk aversion (HARA)
# https://en.wikipedia.org/wiki/Hyperbolic_absolute_risk_aversion
get_upow <- function(certainty, high, low = 0, p = 0.5) {
  log(p) / log(certainty / high)
}

# Example: prefers 20 fish for sure vs. 0/200 gamble at 50-50 odds
get_upow(20, 50)

# upow = 1 (max yield objective)
# upow --> 0 stronger risk aversion and approximates log(catches)

opt_yield <- optim(c(0, 0), get_yield_hcr,
  wt = wt, upow = 1 # max yield
)

# optimize for upow = 0.6 risk-averse (HARA utility)
opt_hara <- optim(c(0, 0), get_yield_hcr,
                  wt = wt, upow = 0.6
)

# extract parameters
lrp_yield <- opt_yield$par[1]
cslope_yield <- opt_yield$par[2]
lrp_hara <- opt_hara$par[1]
cslope_hara <- opt_hara$par[2]

# Create biomass values
xvals <- seq(0.01, sum(ninit * wa * vul), length.out = 200)

# Relative HCRs (your rule)
F_yield <- ifelse(xvals > lrp_yield,
                  cslope_yield * (xvals - lrp_yield) / xvals, 0
)
F_hara <- ifelse(xvals > lrp_hara,
                 cslope_hara * (xvals - lrp_hara) / xvals, 0
)

# linear HCRs
F_yield_lin <- ifelse(xvals > lrp_yield,
                      cslope_yield * (xvals - lrp_yield), 0
)
F_hara_lin <- ifelse(xvals > lrp_hara,
                     cslope_hara * (xvals - lrp_hara), 0
)

# plot
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Panel 1: Relative (space we optimized in)
plot(xvals, F_yield,
     type = "l", lwd = 2, col = "blue",
     xlab = "Vulnerable biomass", ylab = "Fishing mortality (F)",
     main = "Relative HCR: Ft = c*(vB - LRP)/vB", 
     ylim = c(0, max(F_yield, F_hara))
)
lines(xvals, F_hara, col = "red", lwd = 2, lty = 2)
abline(v = lrp_yield, col = "blue", lty = 3)
abline(v = lrp_hara, col = "red", lty = 3)
legend("bottomright",
       legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n"
)

# Panel 2: Linear (space stakeholders and managers want to see)
plot(xvals, F_yield_lin,
     type = "l", lwd = 2, col = "blue",
     xlab = "Vulnerable biomass", ylab = "Fishing mortality (F)",
     main = "Linear HCR: Ft = c*(vB - LRP)", 
     ylim = c(0, max(F_yield_lin, F_hara_lin))
)
lines(xvals, F_hara_lin, col = "red", lwd = 2, lty = 2)
abline(v = lrp_yield, col = "blue", lty = 3)
abline(v = lrp_hara, col = "red", lty = 3)
legend("bottomright",
       legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n"
)

# simulate dynamics
simulate_dynamics <- function(lrp, cslope, wt) {
  log_n <- matrix(-Inf, nrow = n_years, ncol = n_ages)
  log_n[1, ] <- log(ninit)
  ssb <- yield <- vul_bio <- F_out <- numeric(n_years)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)
  for (t in 2:n_years) {
    Ft <- max(0, cslope * (vul_bio[t - 1] - lrp) / vul_bio[t - 1])
    Zt <- Ft * vul + M
    F_out[t] <- Ft
    log_n[t, 1] <- ln_alpha + log(ssb[t - 1]) - br * ssb[t - 1] + wt[t - 1]
    for (a in 2:n_ages) {
      log_n[t, a] <- log_n[t - 1, a - 1] - Zt[a - 1]
    }
    log_n[t, n_ages] <- log(
      exp(log_n[t, n_ages]) +
        exp(log_n[t - 1, n_ages]) * exp(-Zt[n_ages])
    )
    ssb[t] <- sum(exp(log_n[t, ]) * mat * wa)
    vul_bio[t] <- sum(exp(log_n[t, ]) * wa * vul)
    yield[t] <- sum(exp(log_n[t, ]) * wa * Ft * vul / Zt * (1 - exp(-Zt)))
  }
  list(vul_bio = vul_bio, yield = yield, F = F_out)
}

# simulate under each rule
dyn_yield <- simulate_dynamics(lrp_yield, cslope_yield, wt)
dyn_hara  <- simulate_dynamics(lrp_hara, cslope_hara, wt)

# plot biomass, yield, and fishing mortality
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

# 1. Fishing mortality
plot(dyn_yield$F[-1], type = "l", col = "blue", lwd = 2,
     ylab = "Fishing mortality (F)", xlab = "Year", main = "F(t) through time")
lines(dyn_hara$F[-1], col = "red", lwd = 2, lty = 2)
legend("bottomright", legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

# 2. Vulnerable biomass
plot(dyn_yield$vul_bio[-1], type = "l", col = "blue", lwd = 2,
     ylab = "Vulnerable biomass", xlab = "Year", main = "Biomass through time")
lines(dyn_hara$vul_bio[-1], col = "red", lwd = 2, lty = 2)
legend("bottomright", legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

# 3. Yield
plot(dyn_yield$yield[-1], type = "l", col = "blue", lwd = 2,
     ylab = "Yield", xlab = "Year", main = "Yield achieved through time")
lines(dyn_hara$yield[-1], col = "red", lwd = 2, lty = 2)
legend("bottomright", legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

# plot F(t) vs. vulnerable biomass
par(mfrow=c(1,1))
plot(dyn_yield$vul_bio[-1], dyn_yield$F[-1], col = "blue", pch = 16,
     xlab = "Vulnerable biomass", ylab = "Fishing mortality (F)",
     main = "Mortality rate vs. biomass", cex = 0.8)
points(dyn_hara$vul_bio[-1], dyn_hara$F[-1], col = "red", pch = 16, cex = 0.8)
legend("bottomright", legend = c("Max yield", "HARA utility"),
       col = c("blue", "red"), pch = 16, bty = "n")

#----------------------------
# Some stuff worth discussing  
# Dynamic programming solutions from 1975-1990s show these general findings
# Parama, Walters, Hilborn, more Walters, etc.
# Gives you a means to manage a dynamic fishery via a feedback policy
# Which policy is more sustainable? 
# Could you do this with more complex objective functions? yes, but...
# A point about overcapitalization and the so-called precautionary approach