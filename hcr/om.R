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
sdr <- 0.05 # LOOK HERE
set.seed(1)
wt <- rnorm(n_years - 1, 0, sdr)

# objective function
f <- function(par) {
  getAll(data, par)
  Ft <- exp(logF)
  log_n <- matrix(-Inf, n_years, length(mat))
  ssb <- yield <- vul_bio <- numeric(n_years)
  log_n[1, ] <- log(ninit)
  for (t in 1:(n_years - 1)) {
    vul_bio[t] <- sum(exp(log_n[t,]) * wa * vul)
    ssb[t] <- sum(exp(log_n[t,]) * mat * wa)
    Zt <- Ft[t] * vul + M
    yield[t] <- sum(exp(log_n[t,]) * wa * Ft[t] * vul / Zt * (1 - exp(-Zt)))
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
  -mean(yield^upow)
}

# shared data and parameter list
data <- list(
  wt = wt, vul = vul, wa = wa, mat = mat, M = M,
  ln_alpha = ln_alpha, br = br, ninit = ninit,
  upow = NA, n_years = n_years
)
par <- list(logF = rep(log(0.5), n_years - 1))

# yield-maximizing
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
Ft_yield <- exp(opt_yield$par)
rep_yield <- obj_yield$report()

# doone <- function() {
#   nlminb(obj_yield$par + rnorm(length(obj_yield$par),
#     sd = 0.1
#   ), obj_yield$fn, obj_yield$gr)$par
# }
# jit <- replicate(100, doone())
# boxplot(t(jit))

# HARA (risk-averse)
data$upow <- 0.6
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
Ft_hara <- exp(opt_hara$par)
rep_hara <- obj_hara$report()

# setup
vbo <- sum(ninit * wa * vul)
crash_years <- seq(100, n_years, by = 100)
crash_F_indices <- crash_years - 1
keep <- 1:(n_years - 50)
burn_idx <- which((keep - 1) <= 10) # first 10 years as burn in
crash_F_indices <- crash_F_indices[crash_F_indices %in% (keep - 1)]

vb_yield <- rep_yield$vul_bio[keep - 1]
vb_hara <- rep_hara$vul_bio[keep - 1]
Ft_y <- Ft_yield[keep - 1]
Ft_h <- Ft_hara[keep - 1]
highlight_idx <- which((keep - 1) %in% crash_F_indices)

# setup
vbo <- sum(ninit * wa * vul)
crash_years <- seq(100, n_years, by = 100)
keep <- 1:(n_years - 50)

# aligned data
vb_yield <- rep_yield$vul_bio[keep]
vb_hara <- rep_hara$vul_bio[keep]
Ft_y <- Ft_yield[keep]
Ft_h <- Ft_hara[keep]

burn_idx <- which(keep <= 10)
highlight_idx <- which(keep %in% crash_years)

Ut_y <- 1 - exp(-Ft_y)
Ut_h <- 1 - exp(-Ft_h)

# layout: 3 rows of plots
par(mfrow = c(3, 2))
par(mar = c(4, 4.2, 2, 1))

# Panel 1: Relative scale, yield
plot(vb_yield / vbo, Ft_y,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Estimated Ft",
  main = "Yield-maximizing (relative scale)", ylim = c(0, 4.5),
  xlim = c(0, 1)
)
points(vb_yield[highlight_idx] / vbo, Ft_y[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_yield[burn_idx] / vbo, Ft_y[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)

# Panel 2: Relative scale, HARA
plot(vb_hara / vbo, Ft_h,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Estimated Ft",
  main = "Risk-averse (relative scale)", ylim = c(0, 1.15),
  xlim = c(0, 1)
)
points(vb_hara[highlight_idx] / vbo, Ft_h[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_hara[burn_idx] / vbo, Ft_h[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)

# Panel 3: Yield-maximizing (U space)
plot(vb_yield / vbo, Ut_y,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Exploitation rate (U)",
  main = "Yield-maximizing in U space", ylim = c(0, 1), xlim = c(0, 1)
)
points(vb_yield[highlight_idx] / vbo, Ut_y[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_yield[burn_idx] / vbo, Ut_y[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)
abline(h = 0.5, col = "black", lty = 2)

# Panel 4: Risk-averse (U space)
plot(vb_hara / vbo, Ut_h,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Exploitation rate (U)",
  main = "Risk-averse in U space", ylim = c(0, 1), xlim = c(0, 1)
)
points(vb_hara[highlight_idx] / vbo, Ut_h[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_hara[burn_idx] / vbo, Ut_h[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)
abline(h = 0.5, col = "black", lty = 2)

# Panel 5: TAC vs. VB (yield-maximizing)
plot(vb_yield, Ut_y * vb_yield,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Vulnerable biomass", ylab = "Total Allowable Catch (TAC)",
  main = "Yield maximizing policy"
)
points(vb_yield[highlight_idx], Ut_y[highlight_idx] * vb_yield[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_yield[burn_idx], Ut_y[burn_idx] * vb_yield[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)
abline(h = 0.5, col = "black", lty = 2)

# Panel 6: TAC vs. VB (risk-averse)
plot(vb_hara, Ut_h * vb_hara,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Vulnerable biomass", ylab = "Total Allowable Catch (TAC)",
  main = "Risk-averse policy", ylim = c(0, 1), xlim = c(0, 1)
)
points(vb_hara[highlight_idx], Ut_h[highlight_idx] * vb_hara[highlight_idx],
  pch = 4, col = "blue", lwd = 1.2
)
points(vb_hara[burn_idx], Ut_h[burn_idx] * vb_hara[burn_idx],
  pch = 4, col = "red", cex = 0.6
)
abline(v = 1, col = "black", lty = 2)
abline(h = 0.5, col = "black", lty = 2)

# Things worth thinking about
# What do you see, and what do you think is going on?
# What are the points to the far right? What is the OM trying to do?
# How do these patterns change with more recruitment variability?
# Pros and cons of an HCR or OM vs. a basic Fmsy policy?
# What are some ways you could extend this? How does it differ from close-loop
# management strategy evaluation or dynamic programming solutions to the
# feedback policy design problem?
