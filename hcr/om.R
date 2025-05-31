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
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    Zt <- Ft[t - 1] * vul + M
    log_n[t, 1] <- ln_alpha + log(ssb[t - 1]) - br * ssb[t - 1] + wt[t - 1]
    for (a in 2:n_ages) {
      log_n[t, a] <- log_n[t - 1, a - 1] - Zt[a - 1]
    }
    log_n[t, n_ages] <- log(
      exp(log_n[t, n_ages]) + exp(log_n[t - 1, n_ages]) * exp(-Zt[n_ages])
    )
    # modulate population abundance every 100 years or so
    if (t %% 100 == 0) {
      shock <- -10
      log_n[t, ] <- log_n[t, ] + shock # modulate n by e^shock
    }
    n <- exp(log_n[t, ])
    ssb[t] <- sum(n * mat * wa)
    vul_bio[t] <- sum(n * wa * vul)
    yield[t] <- sum(n * wa * Ft[t - 1] * vul / Zt * (1 - exp(-Zt)))
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
keep <- 2:(n_years - 50)
crash_F_indices <- crash_F_indices[crash_F_indices %in% (keep - 1)]

vb_yield <- rep_yield$vul_bio[keep - 1]
vb_hara <- rep_hara$vul_bio[keep - 1]
Ft_y <- Ft_yield[keep - 1]
Ft_h <- Ft_hara[keep - 1]
highlight_idx <- which((keep - 1) %in% crash_F_indices)

layout(matrix(1:4, 2, 2, byrow = TRUE))
par(mar = c(4, 4.2, 2, 1))

# Panel 1: Absolute scale, yield
plot(vb_yield, Ft_y,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Vulnerable biomass", ylab = "Estimated Ft",
  main = "Yield-maximizing (upow = 1)", ylim = c(0, 1.15),
  xlim = c(0, vbo)
)
points(vb_yield[highlight_idx], Ft_y[highlight_idx],
  pch = 4,
  col = "blue", lwd = 1.2
)
abline(v = vbo, col = "black", lty = 2)
legend("topleft",
  legend = "Crash-year Ft", pch = 4,
  col = "blue", bty = "n"
)

# Panel 2: Absolute scale, HARA
plot(vb_hara, Ft_h,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Vulnerable biomass", ylab = "Estimated Ft",
  main = "Risk-averse (upow = 0.6)", ylim = c(0, 1.15),
  xlim = c(0, vbo)
)
points(vb_hara[highlight_idx], Ft_h[highlight_idx],
  pch = 4,
  col = "blue", lwd = 1.2
)
abline(v = vbo, col = "black", lty = 2)
legend("topleft",
  legend = "Crash-year Ft", pch = 4,
  col = "blue", bty = "n"
)

# Panel 3: Relative scale, yield
plot(vb_yield / vbo, Ft_y,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Estimated Ft",
  main = "Yield-maximizing (relative scale)", ylim = c(0, 1.15),
  xlim = c(0, 1)
)
points(vb_yield[highlight_idx] / vbo, Ft_y[highlight_idx],
  pch = 4,
  col = "blue", lwd = 1.2
)
abline(v = 1, col = "black", lty = 2)

# Panel 4: Relative scale, HARA
plot(vb_hara / vbo, Ft_h,
  pch = 19, col = "grey50", cex = 0.6,
  xlab = "Relative vulnerable biomass", ylab = "Estimated Ft",
  main = "Risk-averse (relative scale)", ylim = c(0, 1.15),
  xlim = c(0, 1)
)
points(vb_hara[highlight_idx] / vbo, Ft_h[highlight_idx],
  pch = 4,
  col = "blue", lwd = 1.2
)
abline(v = 1, col = "black", lty = 2)
