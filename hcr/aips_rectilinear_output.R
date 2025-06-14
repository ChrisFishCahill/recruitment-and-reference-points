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

sigma_obs <- 0.1 # observation error SD - often between cv = 0.1-0.3
UMAX <- 0.9 # implementation cap on achievable U
beta <- 80 # softplus approximation

# objective function with precautionary HCR and output control logic
f <- function(par) {
  getAll(data, par)
  lrp <- exp(log_lrp)
  btrigger <- exp(log_btrigger)
  ucap <- 1 / (1 + exp(-logit_ucap)) # logit scale ensures 0 < ucap < 1
  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years - 1)
  Ft <- Ut <- numeric(n_years - 1)
  log_n[1, ] <- log(ninit)
  for (t in 1:(n_years - 1)) {
    # inputs
    ssb[t] <- sum(exp(log_n[t, ]) * mat * wa)
    vul_bio[t] <- vb_true <- sum(exp(log_n[t, ]) * wa * vul)
    # step 0: observe biomass with error
    vb_obs <- vb_true * exp(rnorm(1, 0, sigma_obs))

    # step 1: target U from HCR using observed vb
    slope <- ucap / (btrigger - lrp)
    ramp <- slope * (vb_obs - lrp)
    soft_ramp <- (1 / beta) * log(1 + exp(beta * ramp))
    u_target <- ucap - (1 / beta) * log(1 + exp(beta * (ucap - soft_ramp)))

    # step 2: compute TAC from observed vb and Ftarget
    tac <- vb_obs * exp(-M / 2) * u_target

    # step 3: compute implied U from TAC and true vb
    u_implied <- tac / (vb_true * exp(-M / 2))

    # step 4: apply implementation constraint UMAX via soft cap
    u_realized <- -1 / beta * log(exp(-beta * u_implied) + exp(-beta * UMAX))

    # final Ft and Ut used in dynamics
    Ut[t] <- u_realized
    Ft[t] <- -log(1 - u_realized)
    Zt <- Ft[t] * vul + M
    yield[t] <- sum(exp(log_n[t, ]) * wa * Ft[t] * vul / Zt * (1 - exp(-Zt)))

    # population dynamics
    for (a in 2:n_ages) {
      log_n[t + 1, a] <- log_n[t, a - 1] - Zt[a - 1]
    }
    log_n[t + 1, n_ages] <- log(
      exp(log_n[t + 1, n_ages]) +
        exp(log_n[t, n_ages]) * exp(-Zt[n_ages])
    )
    log_n[t + 1, 1] <- ln_alpha + log(ssb[t]) - br * ssb[t] + wt[t]
    if (t %% 100 == 0) {
      shock <- -10
      log_n[t + 1, ] <- log_n[t + 1, ] + shock
    }
  }
  REPORT(Ut)
  REPORT(yield)
  REPORT(vul_bio)
  REPORT(Ft)
  -mean(yield^upow)
}

data <- list(
  wt = wt, vul = vul, wa = wa, mat = mat, M = M,
  ln_alpha = ln_alpha, br = br, ninit = ninit,
  upow = NA, n_years = n_years
)

# starting values for yield-maximizing HCR
par <- list(
  log_lrp = log(0.2), # lrp
  log_btrigger = log(1.0),
  logit_ucap = qlogis(0.25)
)

# fit yield-maximizing rule
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr)
rep_yield <- obj_yield$report()
logit_ucap <- opt_yield$par["logit_ucap"]
ucap <- 1 / (1 + exp(-logit_ucap))
F <- -log(1 - ucap)
cat(F)

# doone <- function() {
#   nlminb(obj_yield$par + rnorm(length(obj_yield$par),
#                               sd = 0.05
#   ), obj_yield$fn, obj_yield$gr)$par
# }
# jit <- replicate(100, doone())
# boxplot(t(jit))

# fit HARA utility rule
data$upow <- 0.6
par <- list(
  log_lrp = log(0.1),
  log_btrigger = log(1.0),
  logit_ucap = qlogis(0.15)
)
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
rep_hara <- obj_hara$report()
logit_ucap <- opt_hara$par["logit_ucap"]
ucap <- 1 / (1 + exp(-logit_ucap))
F <- -log(1 - ucap)
cat(F)

# doone <- function() {
#   nlminb(obj_hara$par + rnorm(length(obj_hara$par),
#                               sd = 0.05
#   ), obj_hara$fn, obj_hara$gr)$par
# }
# jit <- replicate(100, doone())
# boxplot(t(jit))

# helper ----------------------------------------------------
softplus <- function(z, beta) (1 / beta) * log1p(exp(beta * z))

# extract fitted HCR parameters -----------------------------
## yield fit
lrp_y <- exp(opt_yield$par["log_lrp"])
btrig_y <- exp(opt_yield$par["log_btrigger"])
ucap_y <- 1 / (1 + exp(-opt_yield$par["logit_ucap"]))
slope_y <- ucap_y / (btrig_y - lrp_y)

## HARA fit
lrp_h <- exp(opt_hara$par["log_lrp"])
btrig_h <- exp(opt_hara$par["log_btrigger"])
ucap_h <- 1 / (1 + exp(-opt_hara$par["logit_ucap"]))
slope_h <- ucap_h / (btrig_h - lrp_h)

# biomass sequences for smooth fitted curves ----------------
vb_seq_y <- seq(0.01, max(rep_yield$vul_bio), length.out = 200)
vb_seq_h <- seq(0.01, max(rep_hara$vul_bio), length.out = 200)

# predicted TAC and U(t) using the rectilinear soft-cap rule --
pred_fun <- function(vb, slope, lrp, ucap) {
  ramp <- slope * (vb - lrp)
  soft_ramp <- softplus(ramp, beta)
  Ut <- ucap - softplus(ucap - soft_ramp, beta)
  TAC <- Ut * vb
  list(Ut = Ut, TAC = TAC)
}

pred_y <- pred_fun(vb_seq_y, slope_y, lrp_y, ucap_y)
pred_h <- pred_fun(vb_seq_h, slope_h, lrp_h, ucap_h)

# convenience objects for scatter points --------------------
vb_y <- rep_yield$vul_bio
vb_h <- rep_hara$vul_bio
TAC_y_obs <- rep_yield$yield # TAC = yield here
TAC_h_obs <- rep_hara$yield
Ut_y_obs <- rep_yield$Ut
Ut_h_obs <- rep_hara$Ut

# 2 × 2 panel figure ----------------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# (1) TAC vs vb – yield objective
plot(vb_y, TAC_y_obs,
     pch = 16, col = "dodgerblue3", cex = 0.4,
     xlab = "vulnerable biomass", ylab = "TAC",
     main = "Yield objective", xlim = c(0, 3.0)
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
abline(v = btrig_y, col = "dodgerblue3", lty = 3)

# (2) U(t) vs vb – yield objective
plot(vb_y, Ut_y_obs,
     pch = 16, col = "dodgerblue3", cex = 0.4,
     xlab = "vulnerable biomass", ylab = "U(t)",
     ylim = c(0, 1), xlim = c(0, 3.0),
     main = "Yield objective"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
abline(v = btrig_y, col = "dodgerblue3", lty = 3)
abline(h = UMAX, col = "dodgerblue", lty=3)
lines(vb_seq_y, pred_y$Ut, col = "black", lwd = 2)

# (3) TAC vs vb – HARA objective
plot(vb_h, TAC_h_obs,
     pch = 16, col = "darkorchid", cex = 0.4,
     xlab = "vulnerable biomass", ylab = "TAC",
     ylim = c(0, max(TAC_y_obs)), xlim = c(0, max(TAC_y_obs)),
     main = "HARA objective"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
abline(v = btrig_h, col = "darkorchid", lty = 3)

# (4) U(t) vs vb – HARA objective
plot(vb_h, Ut_h_obs,
     pch = 16, col = "darkorchid", cex = 0.4,
     xlab = "vulnerable biomass", ylab = "U(t)",
     ylim = c(0, 1), xlim = c(0, 3.0),
     main = "HARA objective"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
abline(v = btrig_h, col = "darkorchid", lty = 3)
lines(vb_seq_h, pred_h$Ut, col = "black", lwd = 2)

# plot: time series of Ft, vb, yield
par(mfrow = c(1, 1))
years <- 920:999
t_seq <- years - 919
Ft_y <- obj_yield$report()$"Ft"
Ft_h <- obj_hara$report()$"Ft"
plot(t_seq, Ft_y[years],
     type = "l", col = "dodgerblue3", lwd = 2,
     ylim = c(0, max(Ft_y, Ft_h)), ylab = "Ft", xlab = "Year",
     main = "Ft through time"
)
lines(t_seq, Ft_h[years], col = "darkorchid3", lwd = 2, lty = 2)

plot(t_seq, vb_y[years],
     type = "l", col = "dodgerblue3", lwd = 2,
     ylim = c(0, max(vb_y, vb_h)), ylab = "vb", xlab = "Year",
     main = "Biomass through time"
)
lines(t_seq, vb_h[years], col = "darkorchid3", lwd = 2, lty = 2)
legend("topright",
       legend = c("Yield", "HARA"),
       col = c("dodgerblue3", "darkorchid3"), pch = c(16, 1), lty = 1, bty = "n"
)

plot(t_seq, rep_yield$yield[years],
     type = "l", col = "dodgerblue3", lwd = 2,
     ylab = "Yield", xlab = "Year",
     main = "Yield through time"
)
lines(t_seq, rep_hara$yield[years], col = "darkorchid3", lwd = 2, lty = 2)
legend("topright",
       legend = c("Yield", "HARA"),
       col = c("dodgerblue3", "darkorchid3"), pch = c(16, 1), lty = 1, bty = "n"
)
