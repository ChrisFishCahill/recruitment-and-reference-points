## S^m Ricker fit to Escanaba Lake walleye ---------------------------------
# This parameterization follows Hilborn and Walters 1992
library(FSAdata)
data(WalleyeEL)

S <- WalleyeEL$age5
R <- WalleyeEL$age0

## negative log-likelihood -------------------------------------------------
nll <- function(p, S, R) {
  m <- p[1] # exponent on S
  ln_a <- p[2] # ln(a)
  b <- exp(p[3]) # >0
  sd <- exp(p[4]) # >0
  mu <- m * log(S) + ln_a - b * S
  -sum(dnorm(log(R), mu, sd, log = TRUE))
}

par <- c(
  1, # m start
  log(mean(R / S)), # ln(a)
  log(0.005), # b  (log-scale)
  log(sd(log(R / S))) # sd (log-scale)
)

fit <- optim(par, nll, S = S, R = R, hessian = TRUE)

est <- c(
  m    = fit$par[1],
  ln_a = fit$par[2],
  b    = exp(fit$par[3]),
  sd   = exp(fit$par[4])
)

cat("\nmle estimates:\n")
print(round(est, 4))

## profile m across 0–3 ----------------------------------------------------
mvals <- seq(0, 3, length.out = 100)

nll_prof <- vapply(
  mvals,
  function(mm) {
    optim(
      c(est["ln_a"], log(est["b"]), log(est["sd"])),
      function(p) nll(c(mm, p[1], p[2], p[3]), S, R)
    )$value
  },
  numeric(1)
)

nll_min <- min(nll_prof)
cutoff <- nll_min + 1.92 # 0.5 * χ²₍0.95,1₎

left <- which(mvals < est["m"] & nll_prof <= cutoff)
right <- which(mvals > est["m"] & nll_prof <= cutoff)

m_lo <- if (length(left)) max(mvals[left]) else NA
m_hi <- if (length(right)) min(mvals[right]) else NA

cat("\nprofile 95% ci for m:\n")
print(round(c(lower = m_lo, estimate = est["m"], upper = m_hi), 4))

## fitted values & residuals ----------------------------------------------
xmax <- max(S)
Sg <- seq(1, xmax, length.out = 300)
Rhat <- Sg^est["m"] * exp(est["ln_a"] - est["b"] * Sg)

mu <- est["m"] * log(S) + est["ln_a"] - est["b"] * S
res <- log(R) - mu
fhat <- exp(mu)

## plots -------------------------------------------------------------------
graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

plot(S, R,
  pch = 16, col = "grey40", main = "R vs S",
  xlab = "S", ylab = "R", xlim = c(0, xmax)
)
lines(Sg, Rhat, col = "dodgerblue3", lwd = 2)

plot(S, log(R / S),
  pch = 16, col = "grey40",
  main = "ln(R/S) vs S", xlab = "S", ylab = "ln(R/S)",
  xlim = c(0, xmax)
)
lines(Sg, log(Rhat / Sg), col = "dodgerblue3", lwd = 2)
abline(h = 0, col = "gray70", lty = 2)

plot(fhat, res,
  pch = 16, col = "grey40",
  main = "residuals vs fitted", xlab = "fitted", ylab = "residual"
)
abline(h = 0, col = "gray70", lty = 2)

qqnorm(res, pch = 16, main = "q-q plot")
qqline(res, col = "red")

plot(mvals, nll_prof,
  type = "l", lwd = 2,  
  xlab = "m (fixed)", ylab = "negative log-likelihood",
  main = "profile NLL for m"
)
abline(h = cutoff, col = "gray70", lty = 3)
abline(v = est["m"], col = "gray70", lty = 2)
text(0.5, 0.5, "Escanaba Lake\nwalleye   S^m   Ricker", cex = 1.1)

plot.new()
text(0.5, 0.5, "Escanaba Lake\nwalleye depensatory Ricker", cex = 1.1)
