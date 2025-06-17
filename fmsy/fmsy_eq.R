# define life history vectors
vul <- c(0.05, 0.174, 0.571, 0.899, 1.000, 0.6) # vulnerability at age
wa <- c(0.05, 0.104, 0.145, 0.186, 0.212, 0.253) # weight at age
mat <- c(0.05, 0.20, 0.82, 0.99, 1.00, 1.00) # maturity at age
n_ages <- length(mat) # number of age classes
ages <- 1:n_ages
M <- 0.4 # instantaneous natural mortality

# calculate survivorship in unfished state
lo <- numeric(n_ages)
lo[1] <- 1
for (i in 2:n_ages) lo[i] <- lo[i - 1] * exp(-M)
lo[n_ages] <- lo[n_ages - 1] / (1 - exp(-M)) # plus group

# visualize the life history and vulnerability relationships
par(mfrow = c(2, 2))
plot(vul ~ ages, type = "b", ylab = "vulnerability", xlab = "age")
plot(wa ~ ages, type = "b", ylab = "weight", xlab = "age")
plot(mat ~ ages, type = "b", ylab = "maturity", xlab = "age")
plot(lo ~ ages, type = "b", ylab = "survivorship", xlab = "age")

# set recruitment parameters
Ro <- 10 # unfished recruitment
recK <- 5 # Goodyear compensation ratio
sbro <- sum(lo * mat * wa) # spawning biomass per recruit
ln_alpha <- log(recK / sbro) # log alpha
br <- (ln_alpha + log(sbro)) / (Ro * sbro) # ricker beta

# define objective function to minimize (negative yield)
f <- function(logF) {
  F <- exp(logF)
  n_ages <- length(mat)
  Z <- F * vul + M
  log_lx <- numeric(n_ages)
  log_lx[1] <- log(1) # initialize
  for (a in 2:n_ages) {
    log_lx[a] <- log_lx[a - 1] - Z[a - 1]
  }
  log_lx[n_ages] <- log_lx[n_ages] - log(1 - exp(-Z[n_ages]))
  sbrf <- sum(exp(log_lx) * mat * wa)
  Req <- (ln_alpha + log(sbrf)) / (br * sbrf)
  # equilibrium yield
  YPR <- sum(wa * vul * F / Z * exp(log_lx) * (1 - exp(-Z)))
  Yeq <- YPR * Req
  return(list(neg_yield = -Yeq, sbrf = sbrf, Req = Req))
}

# create a sequence of fishing mortality rates to iterate across
logF_seq <- seq(log(0.001), log(5), length.out = 100)
F_seq <- exp(logF_seq)
yield <- numeric(length(logF_seq))
sbrf <- numeric(length(logF_seq))

for (i in seq_along(logF_seq)) {
  res <- f(logF_seq[i])
  yield[i] <- -res$neg_yield
  sbrf[i] <- res$sbrf
}
Fmsy <- exp(logF_seq[which.max(yield)])
msy <- yield[which.max(yield)]

par(mfrow = c(1, 2))
plot(F_seq, yield, type = "l", xlab = "F", ylab = "Equilibrium yield", lwd = 2)
abline(v = Fmsy, col = "red", lty = 2)
legend("topright", legend = "Fmsy", col = "red", lty = 2, lwd = 2, bty = "n")

plot(F_seq, sbrf / sbro, type = "l", xlab = "F", ylab = "sbrf / sbro", lwd = 2)
abline(v = Fmsy, col = "red", lty = 2)
legend("topright", legend = "Fmsy", col = "red", lty = 2, lwd = 2, bty = "n")

sbrf_at_fmsy <- f(log(Fmsy))$sbrf / sbro
cat("sbrf / sbro at Fmsy:", round(sbrf_at_fmsy, 3), "\n")

#--------------------------------------
# solve via one dimensional optimization
# find the logF that gives the maximum yield
opt <- optimize(
  f = function(logF) f(logF)$neg_yield,
  interval = log(c(0.01, 8))
)
Fmsy <- exp(opt$minimum)
-f(opt$minimum)$neg_yield

#-------------------------------------

# sneaky way to calculate FX%, the F that reduces sbrf to X% of unfished level
X <- 0.35 # X% of sbro

# define a function that returns the relative sbrf / sbro
get_sbrf <- function(logF) f(logF)$sbrf / sbro

# use uniroot to solve for the logF value where sbrf / sbro = target X
FX <- exp(uniroot(function(logF) get_sbrf(logF) - X,
  lower = log(0.001), upper = log(5)
)$root)

par(mfrow = c(1, 2))
plot(F_seq, yield, type = "l", xlab = "F", ylab = "Equilibrium yield", lwd = 2)
abline(v = Fmsy, col = "red", lty = 2)
abline(v = FX, col = "blue", lty = 2)
legend("topright",
  legend = c("Fmsy", "X% sbrf"),
  col = c("red", "blue"), lty = 2, lwd = 2, bty = "n"
)

plot(F_seq, sbrf / sbro,
  type = "l", xlab = "F",
  ylim = c(0.1, 1), ylab = "sbrf / sbro", lwd = 2
)
abline(v = Fmsy, col = "red", lty = 2)
abline(v = FX, col = "blue", lty = 2)
legend("topright",
  legend = c("Fmsy", "X% sbrf"),
  col = c("red", "blue"), lty = 2, lwd = 2, bty = "n"
)

cat("Fmsy:", round(Fmsy, 3), "\n")
cat("sbrf/sbro at Fmsy:", round(sbrf_at_fmsy, 3), "\n")
cat("msy: ", round(-f(opt$minimum)$neg_yield, 3), "\n")
cat("FX%:", round(FX, 3), "\n")

#----------------------------------------------------
# try a bunch of targets
FX_targets <- c(0.2, 0.35, 0.4)
FX_vals <- sapply(FX_targets, function(x) {
  exp(uniroot(function(logF) get_sbrf(logF) - x,
    lower = log(1e-5), upper = log(5)
  )$root)
})
par(mfrow = c(1, 2))
plot(F_seq, yield, type = "l", xlab = "F", ylab = "Equilibrium yield", lwd = 2)
abline(v = FX_vals, col = c("red", "blue", "purple"), lty = 2)
legend("topright",
  legend = paste0("F", FX_targets * 100, "% sbrf / sbro"),
  col = c("red", "blue", "purple"), lty = 2, bty = "n"
)
plot(F_seq, sbrf / sbro,
  type = "l", xlab = "F",
  ylim = c(0.1, 1), ylab = "sbrf / sbro", lwd = 2
)
abline(v = FX_vals, col = c("red", "blue", "purple"), lty = 2)
legend("topright",
  legend = paste0("F", FX_targets * 100, "% sbrf / sbro"),
  col = c("red", "blue", "purple"), lty = 2, bty = "n"
)

#----------------------------------------------------
# maximize yield subject to sbrf/sbro >= 0.65
constrained_f <- function(logF) {
  res <- f(logF)
  if (res$sbrf / sbro < 0.65) {
    return(Inf)
  } # penalize if below threshold
  return(res$neg_yield)
}

opt <- suppressWarnings(optimize(constrained_f, interval = log(c(1e-3, 5))))
F_constrained <- exp(opt$minimum)
F_constrained
