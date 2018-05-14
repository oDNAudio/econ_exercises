library(BMS)

data(datafls)

Y = datafls[, 1]
X = datafls[, 2:ncol(datafls)]

rm(datafls)

ssvs(Y, X)
