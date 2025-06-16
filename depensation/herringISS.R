## depensatory ricker fit to escanaba lake walleye -------------------------
library(FSAdata)
data(HerringISS)
HerringISS <- complete.cases(HeringISS)
S <- HerringISS$ssb[2:50]
R <- HerringISS$rec[2:50]

## negative log-likelihood -------------------------------------------------
nll <- function(p, S, R) {
  a <- p[1]
  b <- exp(p[2])
  h <- exp(p[3])
  sd <- exp(p[4])
  mu <- log(S) + log(S / (h + S)) + a - b * S
  -sum(dnorm(log(R), mu, sd, log = TRUE))
}

par <- c(log(mean(R / S)), log(0.005), log(mean(S) / 2), log(sd(log(R / S))))
fit <- optim(par, nll, S = S, R = R, hessian = TRUE)

## mle estimates -----------------------------------------------------------
est <- c(
  a = fit$par[1],
  b = exp(fit$par[2]),
  h = exp(fit$par[3]),
  sd = exp(fit$par[4])
)

cat("\nmle estimates:\n")
print(round(est, 4))

## profile h across wide range ---------------------------------------------
hvals <- seq(0.1, max(S), length.out = 100)
nll_prof <- numeric(length(hvals))

for (i in seq_along(hvals)) {
  logh <- log(hvals[i])
  nll_prof[i] <- optim(
    c(est["a"], log(est["b"]), log(est["sd"])),
    function(p) nll(c(p[1], p[2], logh, p[3]), S, R)
  )$value
}

## fitted values and residuals ---------------------------------------------
xmax <- max(S)
Sg <- seq(0, xmax, length.out = 300)
Rhat <- (Sg / (est["h"] + Sg)) * Sg * exp(est["a"] - est["b"] * Sg)

mu <- log(S) + log(S / (est["h"] + S)) + est["a"] - est["b"] * S
res <- log(R) - mu
fhat <- exp(mu)

## plots -------------------------------------------------------------------
graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

plot(S, R,
     pch = 16, col = "grey40", main = "R vs S",
     xlab = "S", ylab = "R", xlim = c(0, xmax)
)
lines(Sg, Rhat, col = "dodgerblue3", lwd = 2)
abline(v = est["h"], col = "gray70", lty = 2)

plot(S, log(R / S),
     pch = 16, col = "grey40",
     main = "ln(R/S) vs S", xlab = "S", ylab = "ln(R/S)",
     xlim = c(0, xmax)
)
lines(Sg, log(Rhat / Sg), col = "dodgerblue3", lwd = 2)
abline(v = est["h"], col = "gray70", lty = 2)
abline(h = 0, col = "gray70", lty = 2)

plot(fhat, res,
     pch = 16, col = "grey40",
     main = "residuals vs fitted", xlab = "fitted", ylab = "residual"
)
abline(h = 0, col = "gray70", lty = 2)

qqnorm(res, pch = 16, main = "q-q plot")
qqline(res, col = "red")

plot(hvals, nll_prof,
     type = "l", lwd = 2,
     xlab = "h (fixed)", ylab = "negative log-likelihood",
     main = "profile NLL for h :("
)
abline(v = est["h"], col = "gray70", lty = 2)

plot.new()
text(0.5, 0.5, "Icelandic Herring", cex = 1.1)
