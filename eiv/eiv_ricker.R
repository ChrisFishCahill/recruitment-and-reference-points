set.seed(1)

# --- true model parameters ---
n <- 300
n_sim <- 1000
alpha <- 1.25
beta <- 0.05
log_alpha_true <- log(alpha)

# --- range of true SSB observation error SDs (log-scale) ---
sigma_s_vals <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1.0)

# --- fixed true recruitment observation error SD (log-scale) ---
sigma_r <- 0.5

# --- plot layout ---
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

for (sigma_s in sigma_s_vals) {
  log_alpha_naive <- numeric(n_sim)
  log_alpha_corrected <- numeric(n_sim)

  for (i in 1:n_sim) {
    # -- simulate true values --
    S_true <- runif(n, 1, 10)
    R_true <- alpha * S_true * exp(-beta * S_true)

    # -- simulate observed values with known error SDs --
    S_obs <- S_true * exp(rnorm(n, 0, sigma_s))
    R_obs <- R_true * exp(rnorm(n, 0, sigma_r))

    # -- fit naive model --
    logRS <- log(R_obs / S_obs)
    fit <- lm(logRS ~ S_obs)
    log_alpha_hat <- coef(fit)[1]
    log_alpha_naive[i] <- log_alpha_hat
  }

  # -- plot --
  dens_naive <- density(log_alpha_naive)
  ymax <- max(dens_naive$y)

  plot(dens_naive,
    col = "dodgerblue3", lwd = 2,
    ylim = c(0, ymax), xlim = c(-0.1, 0.9),
    main = bquote(sigma[SSB] == .(sigma_s)),
    xlab = expression(hat(log(alpha)))
  )
  abline(v = log_alpha_true, col = "black", lwd = 2)

  if (sigma_s == sigma_s_vals[1]) {
    legend("topright",
      legend = c("Naive Ricker", "True"), cex = 0.75,
      col = c("dodgerblue3", "black"),
      lwd = 2, bty = "n"
    )
  }
}
