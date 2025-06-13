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

# recruitment params
Ro <- 10
recK <- 5
sbro <- sum(lo * mat * wa)
ln_alpha <- log(recK / sbro)
br <- (ln_alpha + log(sbro)) / (Ro * sbro)
ninit <- Ro * lo

# simulation settings
n_years <- 100
sdr <- 0.6
n_sim <- 1000

# draw wt sequences
set.seed(1)
wt_mat <- matrix(rnorm((n_years - 1) * n_sim, 0, sdr), nrow = n_sim)
# plot 5 example wt sequences
matplot(t(wt_mat[1:5, ]),
  type = "l", lty = 1,
  main = "Examples of Recruitment Deviations (wt)",
  xlab = "Year", ylab = "wt", col = 1:5
)
legend("topright", legend = paste("Sim", 1:5), col = 1:5, lty = 1, bty = "n")

# function to estimate Fmsy from one wt sequence
get_yield <- function(logF, wt) {
  F <- exp(logF)
  Z <- matrix(F * vul + M, nrow = n_years, ncol = n_ages, byrow = TRUE)
  log_n <- matrix(-Inf, nrow = n_years, ncol = n_ages)
  log_n[1, ] <- log(ninit)
  ssb <- yield <- numeric(n_years - 1)
  for (t in 1:(n_years - 1)) {
    ssb[t] <- sum(exp(log_n[t, ]) * mat * wa)
    yield[t] <- sum(exp(log_n[t, ]) * wa * F * vul / Z[t, ] * (1 - exp(-Z[t, ])))
    for (a in 2:n_ages) {
      log_n[t + 1, a] <- log_n[t, a - 1] - Z[t, a - 1]
    }
    log_n[t + 1, n_ages] <- log(
      exp(log_n[t + 1, n_ages]) +
        exp(log_n[t, n_ages]) * exp(-Z[t, n_ages])
    )
    log_n[t + 1, 1] <- ln_alpha + log(ssb[t]) - br * ssb[t] + wt[t]
  }
  -mean(yield)
}

# run simulation
fmsy_vals <- msy_vals <- numeric(n_sim)
for (i in 1:n_sim) {
  wt <- wt_mat[i, ]
  opt <- optimize(function(logF) get_yield(logF, wt),
    interval = log(c(0.01, 5))
  )
  fmsy_vals[i] <- exp(opt$minimum)
  msy_vals[i] <- -get_yield(opt$minimum, wt) # positive MSY
  if (i %% 50 == 0) cat("Completed", i, "of", n_sim, "\n")
}

# plot the stochastic yield curve for one replicate
wt_example <- wt_mat[1, ]
f_seq <- seq(log(0.01), log(5), length.out = 100)
yields <- sapply(f_seq, get_yield, wt = wt_example)

plot(exp(f_seq), -yields,
  type = "l", lwd = 2,
  xlab = "Fishing mortality (F)", ylab = "Yield",
  main = "Stochastic Yield Curve (1 realization)"
)
abline(v = fmsy_vals[1], col = "blue", lty = 2)
legend("topright",
  legend = c("Stochastic Fmsy"),
  col = c("blue"), lty = 2, bty = "n"
)

# plot the same thing but for 3 replicates to visualize
rep_ids <- c(5, 20, 350)
colors <- c("blue", "darkgreen", "orange")

# sequence of F values to evaluate
f_seq <- seq(log(0.01), log(5), length.out = 100)

# compute yields for each replicate
yield_mat <- sapply(rep_ids, function(i) {
  wt_example <- wt_mat[i, ]
  sapply(f_seq, get_yield, wt = wt_example)
})

# plot all yield curves
matplot(exp(f_seq), -yield_mat,
  type = "l", lty = 1, lwd = 2,
  col = colors, xlab = "Fishing mortality (F)", ylab = "Yield",
  main = "Stochastic Yield Curves for 3 Replicates"
)

# add vertical lines for their estimated Fmsy
for (j in seq_along(rep_ids)) {
  abline(v = fmsy_vals[rep_ids[j]], col = colors[j], lty = 2, lwd = 1.5)
}

legend("topright",
  legend = paste("Fmsy for replicate", rep_ids),
  col = colors, lty = 1, lwd = 2, bty = "n"
)

# plot histogram of Fmsy
hist(fmsy_vals,
  breaks = 20,
  main = "Distribution of Fmsy across 1000 simulations",
  xlab = "Estimated Fmsy", col = "lightgray", border = "white"
)

# summarize the median +/- 90% uncertainty interval of Fmsy
quantile(fmsy_vals, c(0.05, 0.5, 0.95))

hist(msy_vals,
  breaks = 40,
  main = "Distribution of MSY across 1000 simulations",
  xlab = "MSY (kg/yr)", col = "lightgray", border = "white"
)

# summarize the median +/- 90% uncertainty interval of Msy
quantile(msy_vals, c(0.05, 0.5, 0.95))

# Class discussion points
# How could we go about exploring the tradeoff between yield and risk
# in this framework? Could we do it in the equilibrium framework?

# Thinkng about recruitment dynamics
# we have done the simplest thing here in terms of wt. What else could be done
# to make the simulation more realistic? How do you think this may change answers?

# Other things that could go wrong, and what might you do to address them?

# Could you adapt this to work with steepness instead of recK?

# Compare the answers you get here with the equilibrium values. What do you notice
# or see?
