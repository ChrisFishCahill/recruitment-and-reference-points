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

  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- numeric(n_years - 1)

  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    log_vb <- log(vul_bio[t - 1])
    delta_log <- log_vb - log(lrp)
    Ft[t - 1] <- cslope * plogis(50 * delta_log) * (1 - exp(-delta_log))

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
      log_n[t, ] <- log_n[t, ] + shock # add shock of e^shock to all ages
    }
    n <- exp(log_n[t, ])
    ssb[t] <- sum(n * mat * wa)
    vul_bio[t] <- sum(n * wa * vul)
    yield[t] <- sum(n * wa * Ft[t - 1] * vul / Zt * (1 - exp(-Zt)))
  }
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

# starting values for smooth HCR parameters
par <- list(log_lrp = log(0.2), log_cslope = log(1))

# fit yield-maximizing rule
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr)
rep_yield <- obj_yield$report()

# doone <- function() {
#   nlminb(obj_yield$par + rnorm(length(obj_yield$par),
#     sd = 0.1
#   ), obj_yield$fn, obj_yield$gr)$par
# }
# jit <- replicate(100, doone())
# boxplot(t(jit))

# fit HARA utility rule
data$upow <- 0.6
par$log_cslope <- log(0.4)
par$log_lrp <- log(1e-10)
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
rep_hara <- obj_hara$report()

# doone <- function() {
#   nlminb(obj_hara$par + rnorm(length(obj_hara$par),
#     sd = 0.1
#   ), obj_hara$fn, obj_hara$gr)$par
# }
# jit <- replicate(100, doone())
# boxplot(t(jit))

# extract fitted parameters
lrp_y <- exp(opt_yield$par["log_lrp"])
cslope_y <- exp(opt_yield$par["log_cslope"])
lrp_h <- exp(opt_hara$par["log_lrp"])
cslope_h <- exp(opt_hara$par["log_cslope"])

# biomass range for drawing the fitted curve
vb_seq <- seq(0.001, max(c(vb_y, vb_h)) * 1.1, length.out = 300)

# fitted Ft curves (smooth log-space HCR)
Ft_y_fit <- cslope_y * plogis(10 * (log(vb_seq) - log(lrp_y))) *
  (1 - exp(-(log(vb_seq) - log(lrp_y))))
Ft_h_fit <- cslope_h * plogis(10 * (log(vb_seq) - log(lrp_h))) *
  (1 - exp(-(log(vb_seq) - log(lrp_h))))

# extract fitted params
lrp_y <- exp(opt_yield$par["log_lrp"])
cslope_y <- exp(opt_yield$par["log_cslope"])
lrp_h <- exp(opt_hara$par["log_lrp"])
cslope_h <- exp(opt_hara$par["log_cslope"])

# vb sequence
vb_y <- rep_yield$vul_bio
vb_h <- rep_hara$vul_bio
Ft_y <- rep_yield$Ft
Ft_h <- rep_hara$Ft
vb_seq <- seq(0.001, max(c(vb_y, vb_h)) * 1.1, length.out = 300)

# Fitted Ft curves
delta_y <- log(vb_seq / lrp_y)
Ft_y_fit <- cslope_y * plogis(10 * delta_y) * (1 - exp(-delta_y))

delta_h <- log(vb_seq / lrp_h)
Ft_h_fit <- cslope_h * plogis(10 * delta_h) * (1 - exp(-delta_h))

# TAC = Ft × vb
TAC_y <- Ft_y_fit * vb_seq
TAC_h <- Ft_h_fit * vb_seq

# 2x2 layout (TAC on top, Ft on bottom)
par(mfrow = c(2, 2), mar = c(4, 4.2, 2, 1))

# 1. TAC vs vb (yield)
plot(vb_seq, TAC_y,
  type = "l", lwd = 2, col = "blue",
  xlab = "Vulnerable biomass", ylab = "Total allowable catch (TAC)",
  main = "TAC = Ft × vb (yield)"
)
abline(v = lrp_y, col = "blue", lty = 3)

# 2. TAC vs vb (HARA)
plot(vb_seq, TAC_h,
  type = "l", lwd = 2, col = "red",
  xlab = "Vulnerable biomass", ylab = "Total allowable catch (TAC)",
  main = "TAC = Ft × vb (HARA)"
)
abline(v = lrp_h, col = "red", lty = 3)

# 3. Ft vs vb (yield)
plot(vb_y[-n_years], Ft_y,
  pch = 16, col = "blue", cex = 0.4,
  xlab = "Vulnerable biomass", ylab = "Fishing mortality Ft",
  main = "Max yield (upow = 1)"
)
abline(v = lrp_y, col = "blue", lty = 3)

# 4. Ft vs vb (HARA)
plot(vb_h[-n_years], Ft_h,
  pch = 16, col = "red", cex = 0.4,
  xlab = "Vulnerable biomass", ylab = "Fishing mortality Ft",
  ylim = c(0, max(Ft_y, Ft_h)),
  main = "HARA utility (upow = 0.6)"
)
abline(v = lrp_h, col = "red", lty = 3)

# ----- plot showing why Ft is nonlinear from linear TAC -----

# choose which parameter set to illustrate (yield or HARA)
lrp_plot <- exp(opt_yield$par["log_lrp"])
cslope_plot <- exp(opt_yield$par["log_cslope"])

# define vb sequence and compute TAC and Ft
vb_seq <- seq(0.01, 2, length.out = 500)
TAC <- cslope_plot * pmax(0, vb_seq - lrp_plot)
Ft_implied <- TAC / vb_seq

# plot
par(mfrow = c(2, 1), mar = c(4, 4.2, 2, 1))

# Panel 1: linear TAC rule
plot(vb_seq, TAC,
  type = "l", lwd = 2, col = "black",
  xlab = "Vulnerable biomass (vb)", ylab = "Total allowable catch (TAC)",
  main = "TAC = c × (vb − LRP) is linear"
)
abline(v = lrp_plot, col = "gray", lty = 3)
text(lrp_plot + 0.1, max(TAC) * 0.8, "LRP", col = "gray")

# Panel 2: resulting nonlinear Ft = TAC / vb
plot(vb_seq, Ft_implied,
  type = "l", lwd = 2, col = "blue",
  xlab = "Vulnerable biomass (vb)", ylab = "Fishing mortality (Ft)",
  main = "Ft = TAC / vb is nonlinear"
)
abline(v = lrp_plot, col = "gray", lty = 3)
text(lrp_plot + 0.1, max(Ft_implied) * 0.8, "LRP", col = "gray")

#- plot dynamics
# select years 920 to 999 (which represent time steps 20 to 100)
years <- 920:999
t_seq <- years - 919 # so it starts at 1 and ends at 80

Ft_y_sub <- rep_yield$Ft[years]
Ft_h_sub <- rep_hara$Ft[years]
vb_y_sub <- rep_yield$vul_bio[years]
vb_h_sub <- rep_hara$vul_bio[years]
yield_y_sub <- rep_yield$yield[years]
yield_h_sub <- rep_hara$yield[years]

# 3-panel plot
par(mfrow = c(3, 1), mar = c(4, 4.2, 2, 1))

# 1. F(t)
plot(t_seq, Ft_y_sub,
  type = "l", lwd = 2, col = "blue",
  ylab = "Fishing mortality (F)", xlab = "Year",
  main = "F(t) through time"
)
lines(t_seq, Ft_h_sub, col = "red", lwd = 2, lty = 2)
legend("topleft",
  legend = c("Max yield", "HARA utility"),
  col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n"
)

# 2. Biomass
plot(t_seq, vb_y_sub,
  type = "l", lwd = 2, col = "blue", ylim = c(0, 3),
  ylab = "Vulnerable biomass", xlab = "Year",
  main = "Biomass through time"
)
lines(t_seq, vb_h_sub, col = "red", lwd = 2, lty = 2)
legend("topleft",
  legend = c("Max yield", "HARA utility"),
  col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n"
)

# 3. Yield
plot(t_seq, yield_y_sub,
  type = "l", lwd = 2, col = "blue",
  ylab = "Yield", xlab = "Year",
  main = "Yield achieved through time"
)
lines(t_seq, yield_h_sub, col = "red", lwd = 2, lty = 2)
legend("topleft",
  legend = c("Max yield", "HARA utility"),
  col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n"
)

# mean yield under yield-maximizing rule
mean_yield_yield <- mean(rep_yield$yield)

# mean yield under HARA utility
mean_yield_hara <- mean(rep_hara$yield)

cat("Mean yield (yield-maximizing):", mean_yield_yield, "\n")
cat("Mean yield (HARA utility):", mean_yield_hara, "\n")

#----------------------------
# Some stuff to discuss
# Dynamic programming solutions from 1975-1990s show these general findings
# Parama, Walters, Hilborn, more Walters, etc.
# Gives you a means to manage a dynamic fishery via a feedback policy
# Which policy is more sustainable?
# Could you do this with more complex objective functions? yes, but likely requires
# simulation (like MSE) or something similar, also...
