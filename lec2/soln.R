# life history vectors
vul <- c(0.05, 0.174, 0.571, 0.899, 1.000, 0.6) # vulnerability at age
vul <- c(0, 0.04, 0.15, 0.1, 0.1, 0.1)
wa <- c(0.05, 0.104, 0.145, 0.186, 0.212, 0.253) # weight at age
mat <- c(0.05, 0.20, 0.82, 0.99, 1.00, 1.00) # maturity at age
M <- 0.4 # instantaneous natural mortality
F <- 0.1 # instantaneous fishing mortality

# visualize the life history and vulnerability relationships
par(mfrow = c(2, 2))
plot(vul, type = "b", ylab = "vulnerability", xlab = "age")
plot(wa, type = "b", ylab = "weight", xlab = "age")
plot(mat, type = "b", ylab = "maturity", xlab = "age")

n_ages <- length(mat) # number of age classes
ages <- 1:n_ages
# set up total instantaneous mortality at age Za vector
Za <- F * vul + M
# calculate survivorship in the fished state
lx <- numeric(n_ages)
lx[1] <- 1 # initialize
for (i in 2:n_ages) lx[i] <- lx[i - 1] * exp(-Za[i - 1]) # calculate lx recursively
lx[n_ages] <- lx[n_ages - 1] / (1 - exp(-Za[n_ages])) # plus group correction
sbrf <- sum(lx * mat * wa) # sumproduct lx, mat, wa
cat("spawner biomass per recruit in the fished condition (sbrf) =", sbrf, "\n")

# calculate survivorship in unfished state (i.e., F = 0)
lo <- numeric(n_ages)
lo[1] <- 1
for (i in 2:n_ages) lo[i] <- lo[i - 1] * exp(-M)
lo[n_ages] <- lo[n_ages - 1] / (1 - exp(-M)) # plus group
sbro <- sum(lo * mat * wa)

# compare lo, lx
par(mfrow = c(1, 1))
plot(lo ~ ages,
  type = "b", pch = 16,
  ylab = "survivorship", col = "dodgerblue3",
  xlab = "age"
)
lines(lx ~ ages, type = "b", pch = 16, col = "darkorchid")
legend("topright",
  legend = c("lo", "lx"),
  col = c("dodgerblue3", "darkorchid"),
  pch = 16,
  lty = 1,
  bty = "n"
)

#-------------------
# Mapping out the equilibrium recruitment function via incidence functions
# set recruitment parameters
Ro <- 10 # unfished recruitment
recK <- 5 # Goodyear compensation ratio
sbro <- sum(lo * mat * wa) # spawning biomass per recruit
alpha <- recK / sbro
br <- (alpha * sbro - 1) / (Ro * sbro) # beverton-holt beta

# sequence of instantaneous fishing mortality rates
Fseq <- seq(0, 8, length.out = 1000)
1 - exp(-Fseq) # finite space exploitation rates (U)
Req <- Seq <- numeric(length(Fseq))
for (i in 1:length(Fseq)) {
  F <- Fseq[i]
  Za <- F * vul + M
  # calculate survivorship in the fished state
  lx <- numeric(n_ages)
  lx[1] <- 1 # initialize
  for (a in 2:n_ages) lx[a] <- lx[a - 1] * exp(-Za[a - 1]) # calculate lx recursively
  lx[n_ages] <- lx[n_ages - 1] / (1 - exp(-Za[n_ages])) # plus group correction
  sbrf <- sum(lx * mat * wa) # sumproduct lx, mat, wa
  Req[i] <- (alpha * sbrf - 1) / (br * sbrf)
  Seq[i] <- Req[i] * sbrf
}
plot(Req ~ Fseq)
abline(lty = 2, h = 0, lwd = 2)
plot(Req ~ Seq)
abline(lty = 2, h = 0, lwd = 2)

#-------------------
# Equilibrium Fmsy/MSY
# set recruitment parameters
Ro <- 10 # unfished recruitment
recK <- 5 # Goodyear compensation ratio
sbro <- sum(lo * mat * wa) # spawning biomass per recruit
alpha <- recK / sbro
br <- (alpha * sbro - 1) / (Ro * sbro) # beverton-holt beta

# sequence of instantaneous fishing mortality rates
Fseq <- seq(0, 8, length.out = 1000)
1 - exp(-Fseq) # finite space exploitation rates (U)
Req <- Seq <- Yeq <- numeric(length(Fseq))
for (i in 1:length(Fseq)) {
  F <- Fseq[i]
  Za <- F * vul + M
  # calculate survivorship in the fished state
  lx <- numeric(n_ages)
  lx[1] <- 1 # initialize
  for (a in 2:n_ages) lx[a] <- lx[a - 1] * exp(-Za[a - 1]) # calculate lx recursively
  lx[n_ages] <- lx[n_ages - 1] / (1 - exp(-Za[n_ages])) # plus group correction
  sbrf <- sum(lx * mat * wa) # sumproduct lx, mat, wa
  Req[i] <- (alpha * sbrf - 1) / (br * sbrf)
  Seq[i] <- Req[i] * sbrf
  YPR <- sum(wa * vul * F / Za * lx * (1 - exp(-Za)))
  Yeq[i] <- YPR * Req[i]
}

par(mfrow=c(1,3))
plot(Req ~ Fseq, main = "R as a function of F")
plot(Req ~ Seq, ylim = c(0, max(Req) + 1), 
     xlim = c(0, max(Seq) + 0.1), main = "Equilibrium R vs. S")
plot(Yeq~Fseq, xlim = c(0,2), ylim = c(0, max(Yeq) + 0.1), 
     main = "Equilibrium yield as a function of F")

# find Fmsy and msy 
idx <- which.max(Yeq)
print(Yeq[idx])
print(Fseq[idx])
