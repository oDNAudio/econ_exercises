# OLS quantities
bols = solve(t(x) %*% x) %*% t(x) %*% y
s2 = t(y - x %*% bols) %*% (y - x %*% bols) / (n - k); s2 = s2[1]
bolscov = s2 * solve(t(x) %*% x)
bolssd = vector("numeric", k)
for(i in 1:k) {
  bolssd[i] = sqrt(bolscov[i, i])
}
v = n - k

# Posterior hyperparameters for Normal-Gamma
xsquare = t(x) %*% x
v1 = v0 + n
capv1inv = capv0inv + xsquare
capv1 = solve(capv1inv)
b1 = capv1 %*% (capv0inv %*% b0 + xsquare %*% bols)
if(det(capv0inv) > 0) {
  v1s12 = v0 %*% s02 + v %*% s2 + t(bols - b0) %*% solve(capv0 + solve(xsquare)) %*% (bols - b0)
} else {
  v1s12 = v0 %*% s02 + v %*% s2
}
v1s12 = v1s12[1]
s12 = v1s12 / v1

bcov = capv1 * v1s12 / (v1 - 2)
bsd = vector("numeric", k)
for(i in 1:k) {
  bsd[i] = sqrt(bcov[i, i])
}

# # Posterior probability of each element of beta being positive + 95% & 99% HPDIs for each
# probpos = vector("numeric", k)
# bhpdi95 = matrix("numeric", k, 2)
# bhpdi99 = matrix("numeric", k, 2)
# 
# invcdf95 = qt(.975, df = v1)
# invcdf99 = qt(.995, df = v1)
# 
# for(i in 1:k) {
#   tnorm = -b1[i] / sqrt(s12 * capv1[i, i])
#   probpos[i] = 1 - pt(tnorm, v1)
#   bhpdi95[i, 1] = b1[i] - invcdf95 * sqrt(s12 * capv1[i, i])
#   bhpdi95[i, 2] = b1[i] + invcdf95 * sqrt(s12 * capv1[i, i])
#   bhpdi99[i, 1] = b1[i] - invcdf99 * sqrt(s12 * capv1[i, i])
#   bhpdi99[i, 2] = b1[i] + invcdf99 * sqrt(s12 * capv1[i, i])
# }

# Posterior mean and variance of error precision
hmean = 1 / s12
hvar = 2 / v1s12
hsd = sqrt(hvar)

# Predictive inference
if(k == 5) {
  xstar = t(c(1, 5000, 2, 2, 1))
  ystarm = xstar %*% b1; ystarm = ystarm[1]
  ystarcapv = (1 + xstar %*% capv1 %*% t(xstar)) * s12; ystarcapv = ystarcapv[1]
  ystarv = ystarcapv * v1 / (v1 - 2)
  ystarsd = sqrt(ystarv)
}

# Log of marginal likelihood if the prior is informative
if(det(capv0inv) != 0) {
  intcon = lgamma(.5 * v1) + .5 * v0 * log(v0 * s02) - lgamma(.5 * v0) - .5 * n * log(pi)
  lmarglik = intcon + .5 * log(det(capv1) / det(capv0)) - .5 * v1 * log(v1s12)
}