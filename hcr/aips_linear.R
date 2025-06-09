library(RTMB)

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

# recruitment
Ro <- 10
recK <- 5
sbro <- sum(lo * mat * wa)
ln_alpha <- log(recK / sbro)
br <- (ln_alpha + log(sbro)) / (Ro * sbro)
ninit <- Ro * lo

# simulation settings
n_years <- 1000
sdr <- 0.6
set.seed(1)
wt <- rnorm(n_years - 1, 0, sdr)

# objective function with smooth HCR
f <- function(par) {
  getAll(data, par)
  lrp <- exp(log_lrp)
  cslope <- exp(log_cslope)
  beta <- 80 # as beta -> inf, approximation -> max(0,z) for softplus approx

  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- Ut <- numeric(n_years - 1)

  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    log_vb <- log(vul_bio[t - 1])
    # linear HCR should specify TAC=max(0,cslope*(vB-lrp)), which implies Ft=-ln(1-TAC/vB)
    # https://en.wikipedia.org/wiki/Softplus Softplus approximation to a
    # non-differentiable function
    Ut[t - 1] <- 1 / beta * log(1 + exp(beta * cslope * (exp(log_vb) - lrp))) / exp(log_vb)
    Ft[t - 1] <- -log(1 - Ut[t - 1])
    Zt <- Ft[t - 1] * vul + M
    log_n[t, 1] <- ln_alpha + log(ssb[t - 1]) - br * ssb[t - 1] + wt[t - 1]
    for (a in 2:n_ages) {
      log_n[t, a] <- log_n[t - 1, a - 1] - Zt[a - 1]
    }
    log_n[t, n_ages] <- log(
      exp(log_n[t, n_ages]) + exp(log_n[t - 1, n_ages]) * exp(-Zt[n_ages])
    )
    if (t %% 100 == 0) {
      shock <- -10
      log_n[t, ] <- log_n[t, ] + shock
    }
    n <- exp(log_n[t, ])
    ssb[t] <- sum(n * mat * wa)
    vul_bio[t] <- sum(n * wa * vul)
    yield[t] <- sum(n * wa * Ft[t - 1] * vul / Zt * (1 - exp(-Zt)))
  }
  REPORT(yield)
  REPORT(vul_bio)
  REPORT(Ft)
  REPORT(Ut)
  -mean(yield^upow)
}

data <- list(
  wt = wt, vul = vul, wa = wa, mat = mat, M = M,
  ln_alpha = ln_alpha, br = br, ninit = ninit,
  upow = NA, n_years = n_years
)

# starting values for smooth HCR parameters
par <- list(log_lrp = log(0.2), log_cslope = log(1))

# fit yield-maximizing rule
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr)
rep_yield <- obj_yield$report()

# fit HARA utility rule
data$upow <- 0.6
par$log_cslope <- log(0.3)
par$log_lrp <- log(1.5)
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
rep_hara <- obj_hara$report()

# ----- helper function -----
softplus <- function(z, beta) (1 / beta) * log1p(exp(beta * z))
beta <- 80 # same as in model

# fitted parameters
lrp_y <- exp(opt_yield$par["log_lrp"])
cslope_y <- exp(opt_yield$par["log_cslope"])
lrp_h <- exp(opt_hara$par["log_lrp"])
cslope_h <- exp(opt_hara$par["log_cslope"])

# simulation outputs
vb_y <- rep_yield$vul_bio
vb_h <- rep_hara$vul_bio
Ft_y <- rep_yield$Ft
Ft_h <- rep_hara$Ft

# Ut and TAC (simulation)
Ut_y <- 1 - exp(-Ft_y)
Ut_h <- 1 - exp(-Ft_h)
TAC_y_obs <- Ut_y * vb_y[-n_years]
TAC_h_obs <- Ut_h * vb_h[-n_years]

# Fitted HCR predictions
TAC_y_pred <- softplus(cslope_y * (vb_y[-n_years] - lrp_y), beta)
Ft_y_pred <- TAC_y_pred / vb_y[-n_years]
Ut_y_pred <- 1 - exp(-Ft_y_pred)

TAC_h_pred <- softplus(cslope_h * (vb_h[-n_years] - lrp_h), beta)
Ft_h_pred <- TAC_h_pred / vb_h[-n_years]
Ut_h_pred <- 1 - exp(-Ft_h_pred)

# sort vb and predictions before drawing lines
ix_y <- order(vb_y[-n_years])
ix_h <- order(vb_h[-n_years])

# plot layout
par(mfrow = c(2, 2), mar = c(4, 4.2, 2, 1))

# 1. TAC vs vb (Yield)
plot(vb_y[-1], TAC_y_obs,
     type = "p", pch = 16, col = "dodgerblue3", cex = 0.4,
     xlab = "Vulnerable biomass (vb)",
     ylab = "TAC = U × vb",
     main = "TAC vs vb (Yield)"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_y[-n_years][ix_y], TAC_y_pred[ix_y], col = "dodgerblue3", lwd = 2)

# 2. TAC vs vb (HARA)
plot(vb_h[-1], TAC_h_obs,
     type = "p", pch = 16, col = "darkorchid3", cex = 0.4,
     xlab = "Vulnerable biomass (vb)",
     ylab = "TAC = U × vb",
     main = "TAC vs vb (HARA)"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_h[-n_years][ix_h], TAC_h_pred[ix_h], col = "darkorchid3", lwd = 2)

# clean input sequence for vb
vb_seq_y <- seq(min(vb_y), max(vb_y), length.out = 200)
vb_seq_h <- seq(min(vb_h), max(vb_h), length.out = 200)

Ft_y_smooth <- -log(1 - softplus(cslope_y * (vb_seq_y - lrp_y), beta) / vb_seq_y)
Ft_h_smooth <- -log(1 - softplus(cslope_h * (vb_seq_h - lrp_h), beta) / vb_seq_h)

# 3. Ft vs vb (Yield)
plot(vb_y[-1], Ft_y,
     type = "p", pch = 16, col = "dodgerblue3", cex = 0.4,
     xlab = "Vulnerable biomass (vb)",
     ylab = "Fishing mortality (Ft)",
     main = "Ft vs vb (Yield)"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_seq_y, Ft_y_smooth, col = "dodgerblue3", lwd = 2)

# 4. Ft vs vb (HARA)
plot(vb_h[-1], Ft_h,
     type = "p", pch = 16, col = "darkorchid3", cex = 0.4,
     xlab = "Vulnerable biomass (vb)",
     ylab = "Fishing mortality (Ft)",
     main = "Ft vs vb (HARA)",
     ylim = c(0, max(Ft_y, Ft_h))
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_seq_h, Ft_h_smooth, col = "darkorchid3", lwd = 2)

#----------------------------
# Some stuff to discuss
# Dynamic programming solutions from 1975-1990s show these general findings
# Parama, Walters, Hilborn, more Walters, etc.
# Gives you a means to manage a dynamic fishery via a feedback policy
# Which policy is more sustainable?
# Could you do this with more complex objective functions? Yes, but may require
# simulation (like MSE) or something similar, also...
