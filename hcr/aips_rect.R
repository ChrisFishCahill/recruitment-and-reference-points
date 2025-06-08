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
  csl <- exp(log_csl)
  bhalf <- exp(log_bhalf)
  fmax <- exp(log_fmax)

  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- numeric(n_years - 1)

  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)

  for (t in 2:n_years) {
    log_vb <- log(vul_bio[t - 1])
    Ft[t - 1] <- fmax / (1 + exp(-csl * (exp(log_vb) - bhalf)))

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
  -mean(yield^upow)
}

data <- list(
  wt = wt, vul = vul, wa = wa, mat = mat, M = M,
  ln_alpha = ln_alpha, br = br, ninit = ninit,
  upow = NA, n_years = n_years
)

# starting values for smooth HCR parameters
par <- list(log_csl = log(0.2), log_bhalf = log(5), log_fmax = log(0.8))

# fit yield-maximizing rule
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr)
rep_yield <- obj_yield$report()

doone <- function() {
  nlminb(obj_yield$par + rnorm(length(obj_yield$par),
    sd = 0.1
  ), obj_yield$fn, obj_yield$gr)$par
}
jit <- replicate(100, doone())
boxplot(t(jit))

plot(obj_yield$report()$`Ft` ~ obj_yield$report()$`vul_bio`[-1])

# fit HARA utility rule
data$upow <- 0.6
par <- list(log_csl = log(1), log_bhalf = log(1), log_fmax = log(1))
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
  control = list(eval.max = 10000, iter.max = 10000)
)
rep_hara <- obj_hara$report()

doone <- function() {
  nlminb(obj_hara$par + rnorm(length(obj_hara$par),
    sd = 0.1
  ), obj_hara$fn, obj_hara$gr)$par
}
jit <- replicate(100, doone())
boxplot(t(jit))
plot(obj_hara$report()$`Ft` ~ obj_hara$report()$`vul_bio`[-1])
