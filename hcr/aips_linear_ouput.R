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

sigma_obs <- 0.2 # observation error SD
UMAX <- 0.9 # implementation cap on achievable U
beta <- 80 # softplus approximation

# objective function with smooth HCR
f <- function(par) {
  getAll(data, par)
  lrp <- exp(log_lrp)
  cslope <- exp(log_cslope)
  log_n <- matrix(-Inf, n_years, length(mat))
  ssb <- yield <- vul_bio <- numeric(n_years - 1)
  Ft <- Ut <- numeric(n_years - 1)
  log_n[1, ] <- log(ninit)
  for (t in 1:(n_years - 1)) {
    vul_bio[t] <- vb_true <- sum(exp(log_n[t, ]) * wa * vul)
    ssb[t] <- sum(exp(log_n[t, ]) * mat * wa)
    # step 0: observe biomass with error
    vb_obs <- vb_true * exp(rnorm(1, 0, sigma_obs)) # observation error

    # step 1: compute target TAC using softened linear ramp
    tac_target <- (1 / beta) * log1p(exp(beta * cslope * (vb_obs - lrp))) # softplus on VB_obs

    # step 2: derive implied exploitation rate based on true biomass
    u_implied <- tac_target / (vb_true * exp(-M / 2)) # output control → backward to U

    # step 3: apply soft implementation constraint (Umax)
    u_realized <- -1 / beta * log(exp(-beta * u_implied) + exp(-beta * UMAX)) # soft min

    # step 4: convert to F and use in population dynamics
    Ut[t] <- u_realized
    Ft[t] <- -log(1 - u_realized)
    Zt <- Ft[t] * vul + M

    yield[t] <- sum(exp(log_n[t, ]) * wa * Ft[t] * vul / Zt * (1 - exp(-Zt)))
    ## survival & ageing
    for (a in 2:n_ages) {
      log_n[t + 1, a] <- log_n[t, a - 1] - Zt[a - 1]
    }
    ## plus-group
    log_n[t + 1, n_ages] <- log(
      exp(log_n[t + 1, n_ages]) +
        exp(log_n[t, n_ages]) * exp(-Zt[n_ages])
    )
    ## recruitment
    log_n[t + 1, 1] <- ln_alpha + log(ssb[t]) - br * ssb[t] + wt[t]
    ## occasional shock of e^-shock
    if (t %% 100 == 0) {
      shock <- -10
      log_n[t + 1, ] <- log_n[t + 1, ] + shock
    }
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

# fitted parameters
lrp_y <- exp(opt_yield$par["log_lrp"])
cslope_y <- exp(opt_yield$par["log_cslope"])
lrp_h <- exp(opt_hara$par["log_lrp"])
cslope_h <- exp(opt_hara$par["log_cslope"])

# plots ---------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# softplus and fitted curves
softplus <- function(z, beta) (1 / beta) * log1p(exp(beta * z))
vb_grid <- seq(0, 3.2, length.out = 200)
tac_y <- softplus(cslope_y * (vb_grid - lrp_y), beta)
tac_h <- softplus(cslope_h * (vb_grid - lrp_h), beta)

# original plots
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(rep_yield$vul_bio, rep_yield$yield,
  pch = 16, col = "dodgerblue3", cex = 0.4,
  xlab = "vulnerable biomass", ylab = "TAC",
  main = "Yield objective"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_grid, pmax(0, tac_y), col = "black", lwd = 2)

plot(rep_yield$vul_bio, rep_yield$yield / rep_yield$vul_bio,
  pch = 16, col = "dodgerblue3", cex = 0.4,
  xlab = "vulnerable biomass", ylab = "U(t)",
  ylim = c(0, 1), main = "Yield objective"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_grid, pmax(0, tac_y) / vb_grid, col = "black", lwd = 2)

plot(rep_hara$vul_bio, rep_hara$yield,
  pch = 16, col = "darkorchid", cex = 0.4,
  xlab = "vulnerable biomass", ylab = "TAC",
  ylim = c(0, 3), xlim = c(0, 3.0), main = "HARA objective"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_grid, pmax(0, tac_h), col = "black", lwd = 2)

plot(rep_hara$vul_bio, rep_hara$yield / rep_hara$vul_bio,
  pch = 16, col = "darkorchid", cex = 0.4,
  xlab = "vulnerable biomass", ylab = "U(t)",
  ylim = c(0, 1), main = "HARA objective"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_grid, pmax(0, tac_h) / vb_grid, col = "black", lwd = 2)
