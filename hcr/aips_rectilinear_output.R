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

sigma_obs <- 1e-3 # observation error SD
UMAX <- 0.9 # implementation cap on achievable U
beta <- 80 # softplus approximation

# objective function with precautionary HCR and output control logic
f <- function(par) {
  getAll(data, par)
  lrp <- exp(log_lrp)
  btrigger <- exp(log_btrigger)
  ucap <- 1 / (1 + exp(-logit_ucap)) # HCR cap on target U

  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- Ut <- numeric(n_years - 1)

  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    vb_true <- vul_bio[t - 1]

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
    Ut[t - 1] <- u_realized
    Ft[t - 1] <- -log(1 - u_realized)

    # population dynamics
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

# ---- helper ----
softplus <- function(z, beta) (1 / beta) * log1p(exp(beta * z))
beta <- 80

# extract fitted parameters
lrp_y <- exp(opt_yield$par["log_lrp"])
btrig_y <- exp(opt_yield$par["log_btrigger"])
ucap_y <- 1 / (1 + exp(-opt_yield$par["logit_ucap"]))
slope_y <- ucap_y / (btrig_y - lrp_y)

lrp_h <- exp(opt_hara$par["log_lrp"])
btrig_h <- exp(opt_hara$par["log_btrigger"])
ucap_h <- 1 / (1 + exp(-opt_hara$par["logit_ucap"]))
slope_h <- ucap_h / (btrig_h - lrp_h)

# biomass vectors
vb_y <- rep_yield$vul_bio
vb_h <- rep_hara$vul_bio
vb_seq <- seq(0.001, max(c(vb_y, vb_h)) * 1.05, length.out = 300)

# smooth predicted Ut
ramp_y <- slope_y * (vb_seq - lrp_y)
soft_ramp_y <- softplus(ramp_y, beta)
Ut_y_pred <- ucap_y - softplus(ucap_y - soft_ramp_y, beta)

ramp_h <- slope_h * (vb_seq - lrp_h)
soft_ramp_h <- softplus(ramp_h, beta)
Ut_h_pred <- ucap_h - softplus(ucap_h - soft_ramp_h, beta)

# align Ut and vb lengths
n_y <- min(length(vb_y) - 1, length(rep_yield$Ut))
n_h <- min(length(vb_h) - 1, length(rep_hara$Ut))

# plot: Ut vs vb
par(mfrow = c(2, 2), mar = c(4, 4.2, 2, 1))
plot(vb_y[-1][1:n_y], rep_yield$Ut[1:n_y],
  col = "dodgerblue3", pch = 16, cex = 0.4,
  xlab = "Vulnerable biomass", ylab = "Exploitation rate (Ut)",
  main = "Ut vs vb"
)
points(vb_h[-1][1:n_h], rep_hara$Ut[1:n_h], col = "darkorchid3", pch = 1, cex = 0.4)
lines(vb_seq, Ut_y_pred, col = "dodgerblue3", lwd = 2)
lines(vb_seq, Ut_h_pred, col = "darkorchid3", lwd = 2)
legend("bottomright",
  legend = c("Yield", "HARA"),
  col = c("dodgerblue3", "darkorchid3"), pch = c(16, 1), lty = 1, bty = "n"
)

# plot: Ft vs vb
Ft_y <- rep_yield$Ft[1:n_y]
Ft_h <- rep_hara$Ft[1:n_h]
Ft_y_pred <- -log(1 - Ut_y_pred)
Ft_h_pred <- -log(1 - Ut_h_pred)

plot(vb_y[-1][1:n_y], Ft_y,
  col = "dodgerblue3", pch = 16, cex = 0.4,
  xlab = "Vulnerable biomass", ylab = "Fishing mortality (Ft)",
  main = "Ft vs vb"
)
points(vb_h[-1][1:n_h], Ft_h, col = "darkorchid3", pch = 1, cex = 0.4)
lines(vb_seq, Ft_y_pred, col = "dodgerblue3", lwd = 2)
lines(vb_seq, Ft_h_pred, col = "darkorchid3", lwd = 2)

# plot: TAC vs vb
TAC_y <- rep_yield$Ut[1:n_y] * vb_y[-1][1:n_y]
TAC_h <- rep_hara$Ut[1:n_h] * vb_h[-1][1:n_h]
TAC_y_pred <- Ut_y_pred * vb_seq
TAC_h_pred <- Ut_h_pred * vb_seq

plot(vb_y[-1][1:n_y], TAC_y,
  col = "dodgerblue3", pch = 16, cex = 0.4,
  xlab = "Vulnerable biomass", ylab = "TAC",
  main = "TAC vs vb"
)
points(vb_h[-1][1:n_h], TAC_h, col = "darkorchid3", pch = 1, cex = 0.4)
lines(vb_seq, TAC_y_pred, col = "dodgerblue3", lwd = 2)
lines(vb_seq, TAC_h_pred, col = "darkorchid3", lwd = 2)

# plot: time series of Ft
years <- 920:999
t_seq <- years - 919
plot(t_seq, rep_yield$Ft[years],
  type = "l", col = "dodgerblue3", lwd = 2,
  ylim = c(0, max(rep_yield$Ft, rep_hara$Ft)), ylab = "Ft", xlab = "Year",
  main = "Ft through time"
)
lines(t_seq, rep_hara$Ft[years], col = "darkorchid3", lwd = 2, lty = 2)

# biomass and yield through time
par(mfrow = c(1, 1))
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
  ylim = c(0, max(rep_yield$yield, rep_hara$yield)),
  ylab = "Yield", xlab = "Year",
  main = "Yield through time"
)
lines(t_seq, rep_hara$yield[years], col = "darkorchid3", lwd = 2, lty = 2)
legend("topright",
  legend = c("Yield", "HARA"),
  col = c("dodgerblue3", "darkorchid3"), pch = c(16, 1), lty = 1, bty = "n"
)
