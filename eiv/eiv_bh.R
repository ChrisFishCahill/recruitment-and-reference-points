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
  log_alpha_nls <- numeric(n_sim)

  for (i in 1:n_sim) {
    # -- simulate true values --
    S_true <- runif(n, 1, 10)
    R_true <- (alpha * S_true) / (1 + beta * S_true)

    # -- simulate observed values with lognormal error --
    S_obs <- S_true * exp(rnorm(n, 0, sigma_s))
    R_obs <- R_true * exp(rnorm(n, 0, sigma_r))

    # -- fit Beverton-Holt model via nls --
    fit <- try(
      nls(R_obs ~ alpha_hat * S_obs / (1 + beta_hat * S_obs),
        start = list(alpha_hat = 1, beta_hat = 0.01),
        control = list(warnOnly = TRUE)
      ),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {
      log_alpha_hat <- log(coef(fit)["alpha_hat"])
      # apply Jensen's inequality correction
      log_alpha_nls[i] <- log_alpha_hat - (sigma_r^2) / 2
    } else {
      log_alpha_nls[i] <- NA
    }
  }

  # -- plot --
  dens_nls <- density(na.omit(log_alpha_nls))
  ymax <- max(dens_nls$y)

  plot(dens_nls,
    col = "dodgerblue3", lwd = 2,
    ylim = c(0, ymax), xlim = c(-0.1, 4),
    main = bquote(sigma[SSB] == .(sigma_s)),
    xlab = expression(hat(log(alpha)))
  )
  abline(v = log_alpha_true, col = "black", lwd = 2)

  if (sigma_s == sigma_s_vals[1]) {
    legend("topright",
      legend = c("Naive Beverton-Holt", "True"), cex = 0.75,
      col = c("dodgerblue3", "black"),
      lwd = 2, bty = "n"
    )
  }
}
