library(BMS)

source("exercise_3_ssvs.R")

data(datafls)

Y = datafls[, 1]
X = datafls[, 2:ncol(datafls)]
rm(datafls)

tau0s = c(1, 1e-2, 1e-5, 1e-15)
tau1_scales = c(100, 1000, 10000)

pips = vector("list", length(tau0s) * length(tau1_scales))

i = 1
for(scale in tau1_scales) {
  for(tau in tau0s) {
    pips[[i]] = ssvs(Y, X, tau0 = tau, tau1 = tau * scale)[[1]]
    names(pips)[i] = paste0(tau, " & *", scale)
    i = i + 1
  }
}

pips_z = vector("list", length(tau0s) * length(tau1_scales))

i = 1
for(scale in tau1_scales) {
  for(tau in tau0s) {
    pips_z[[i]] = ssvs(Y, X, tau0 = tau, tau1 = tau * scale, standardise = TRUE)[[1]]
    names(pips_z)[i] = paste0(tau, " & *", scale)
    i = i + 1
  }
}

0.00001 + 1e-1
