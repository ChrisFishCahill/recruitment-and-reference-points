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

  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- Ut <- numeric(n_years - 1)

  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    vb_true <- vul_bio[t - 1]
    # step 0: observe biomass with error
    vb_obs <- vb_true * exp(rnorm(1, 0, sigma_obs)) # observation error

    # step 1: compute target TAC using softened linear ramp
    log_vb_obs <- log(vb_obs)
    tac_target <- (1 / beta) * log1p(exp(beta * cslope * (vb_obs - lrp))) # softplus on VB_obs

    # step 2: derive implied exploitation rate based on true biomass
    u_implied <- tac_target / (vb_true * exp(-M / 2)) # output control → backward to U

    # step 3: apply soft implementation constraint (Umax)
    u_realized <- -1 / beta * log(exp(-beta * u_implied) + exp(-beta * UMAX)) # soft min

    # step 4: convert to F and use in population dynamics
    Ut[t - 1] <- u_realized
    Ft[t - 1] <- -log(1 - u_realized)
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
## ----- helper (soft-min with UMAX) -----------------
soft_min <- function(x, umax, beta) {
  umax - (1 / beta) * log1p(exp(beta * (umax - x)))
}

## ----- PREPARE VECTORS ALIGNED WITH  F_{t-1}  ------
vb_y_plot <- vb_y[-n_years] # drop last year → length = n_years-1
vb_h_plot <- vb_h[-n_years]

Ft_y <- rep_yield$Ft # already length n_years-1
Ft_h <- rep_hara$Ft

## ----- PREDICTED CURVES  (Yield) -------------------
tac_y_hat <- softplus(cslope_y * (vb_y_plot - lrp_y), beta) # TAC*
u_y_hat <- tac_y_hat / (vb_y_plot * exp(-M / 2)) # implied U
u_y_hat <- soft_min(u_y_hat, UMAX, beta) # cap at UMAX
Ft_y_hat <- -log(1 - u_y_hat) # to F

## ----- PREDICTED CURVES  (HARA) --------------------
tac_h_hat <- softplus(cslope_h * (vb_h_plot - lrp_h), beta)
u_h_hat <- tac_h_hat / (vb_h_plot * exp(-M / 2))
u_h_hat <- soft_min(u_h_hat, UMAX, beta)
Ft_h_hat <- -log(1 - u_h_hat)

## ----- CLEAN x-sequences FOR SMOOTH LINES ----------
vb_seq_y <- seq(min(vb_y_plot), max(vb_y_plot), length.out = 200)
vb_seq_h <- seq(min(vb_h_plot), max(vb_h_plot), length.out = 200)

Ft_y_smooth <- {
  tac <- softplus(cslope_y * (vb_seq_y - lrp_y), beta)
  u <- tac / (vb_seq_y * exp(-M / 2))
  u <- soft_min(u, UMAX, beta)
  -log(1 - u)
}

Ft_h_smooth <- {
  tac <- softplus(cslope_h * (vb_seq_h - lrp_h), beta)
  u <- tac / (vb_seq_h * exp(-M / 2))
  u <- soft_min(u, UMAX, beta)
  -log(1 - u)
}

## ----- PLOT LAYOUT ---------------------------------
par(mfrow = c(2, 2), mar = c(4, 4.2, 2, 1))

# 1. TAC vs vb (Yield)
plot(vb_y_plot, tac_y_hat,
  pch = 16, col = "dodgerblue3", cex = 0.4,
  xlab = "Vulnerable biomass (VB)",
  ylab = "TAC",
  main = "TAC vs VB (Yield)"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_seq_y, softplus(cslope_y * (vb_seq_y - lrp_y), beta),
  col = "dodgerblue3", lwd = 2
)

# 2. TAC vs vb (HARA)
plot(vb_h_plot, tac_h_hat,
  pch = 16, col = "darkorchid3", cex = 0.4,
  xlab = "Vulnerable biomass (VB)",
  ylab = "TAC",
  main = "TAC vs VB (HARA)"
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_seq_h, softplus(cslope_h * (vb_seq_h - lrp_h), beta),
  col = "darkorchid3", lwd = 2
)

# 3. Ft vs vb (Yield)
plot(vb_y_plot, Ft_y,
  pch = 16, col = "dodgerblue3", cex = 0.4,
  xlab = "Vulnerable biomass (VB)",
  ylab = "Fishing mortality (F)",
  main = "F vs VB (Yield)"
)
abline(v = lrp_y, col = "dodgerblue3", lty = 3)
lines(vb_seq_y, Ft_y_smooth, col = "dodgerblue3", lwd = 2)

# 4. Ft vs vb (HARA)
plot(vb_h_plot, Ft_h,
  pch = 16, col = "darkorchid3", cex = 0.4,
  xlab = "Vulnerable biomass (VB)",
  ylab = "Fishing mortality (F)",
  main = "F vs VB (HARA)",
  ylim = c(0, max(Ft_y, Ft_h, Ft_y_smooth, Ft_h_smooth))
)
abline(v = lrp_h, col = "darkorchid3", lty = 3)
lines(vb_seq_h, Ft_h_smooth, col = "darkorchid3", lwd = 2)
