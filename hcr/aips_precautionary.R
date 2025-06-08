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

# objective function with precautionary HCR
f <- function(par) {
  getAll(data, par)
  blim <- exp(log_blim)
  btrigger <- exp(log_btrigger)
  ucap <- 1 / (1 + exp(-logit_ucap))  # logit scale ensures 0 < ucap < 1
  beta <- 80
  
  log_n <- matrix(-Inf, n_years, n_ages)
  ssb <- yield <- vul_bio <- numeric(n_years)
  Ft <- Ut <- numeric(n_years - 1)
  
  log_n[1, ] <- log(ninit)
  ssb[1] <- sum(exp(log_n[1, ]) * mat * wa)
  vul_bio[1] <- sum(exp(log_n[1, ]) * wa * vul)
  
  for (t in 2:n_years) {
    vb <- vul_bio[t - 1]
    slope <- ucap / (btrigger - blim)
    ramp <- slope * (vb - blim)
    soft_ramp <- (1 / beta) * log(1 + exp(beta * ramp))
    Ut[t] <- ucap - (1 / beta) * log(1 + exp(beta * (ucap - soft_ramp)))
    Ft[t - 1] <- -log(1 - Ut[t])
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
  REPORT(Ut)
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

# starting values for yield-maximizing HCR
par <- list(
  log_blim = log(0.2),
  log_btrigger = log(1.0),
  logit_ucap = qlogis(0.25)  # 25% exploitation cap
)

# fit yield-maximizing rule
data$upow <- 1
obj_yield <- MakeADFun(f, par, data = data)
opt_yield <- nlminb(obj_yield$par, obj_yield$fn, obj_yield$gr)
rep_yield <- obj_yield$report()
logit_ucap <- opt_yield$par["logit_ucap"]
ucap <- 1 / (1 + exp(-logit_ucap))
F <- -log(1 - ucap)
cat(F)

doone <- function() {
  nlminb(obj_yield$par + rnorm(length(obj_yield$par),
                              sd = 0.05
  ), obj_yield$fn, obj_yield$gr)$par
}
jit <- replicate(100, doone())
boxplot(t(jit))

# fit HARA utility rule
data$upow <- 0.6
par <- list(
  log_blim = log(0.3),
  log_btrigger = log(1.0),
  logit_ucap = qlogis(0.15)  # 15% exploitation cap for HARA
)
obj_hara <- MakeADFun(f, par, data = data)
opt_hara <- nlminb(obj_hara$par, obj_hara$fn, obj_hara$gr,
                   control = list(eval.max = 10000, iter.max = 10000))
rep_hara <- obj_hara$report()
logit_ucap <- opt_hara$par["logit_ucap"]
ucap <- 1 / (1 + exp(-logit_ucap))
F <- -log(1 - ucap)
cat(F)

doone <- function() {
  nlminb(obj_hara$par + rnorm(length(obj_hara$par),
                              sd = 0.05
  ), obj_hara$fn, obj_hara$gr)$par
}
jit <- replicate(100, doone())
boxplot(t(jit))

# comparison plot
#pdf("compare_hcr_Ut_vs_vulbio.pdf", width = 7, height = 5)
plot(rep_yield$vul_bio, rep_yield$Ut, pch = 16, col = "blue",
     xlab = "Vulnerable biomass", ylab = "Exploitation rate (U)",
     main = "HCR Comparison: Yield vs HARA")
points(rep_hara$vul_bio, rep_hara$Ut, pch = 1, col = "darkred")
legend("bottomright", legend = c("Yield-optimizing", "HARA-optimizing"),
       col = c("blue", "darkred"), pch = c(16, 1), bty = "n")
#dev.off()

