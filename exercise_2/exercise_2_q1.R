library(readr)

setwd("Aufgaben")
data <- read_table2("exercise_2_data.txt")

n = nrow(data)
y = as.matrix(data[, 1])
x = as.matrix(data.frame(rep(1, n), data[, 2:5]))
k = 5

v0 = 5
b0 = c(0, 10, 5000, 10000, 10000)
s02 = 1 / 4.0e-8
capv0 = 2.4 * diag(k)
capv0[2, 2] = 6e-7
capv0[3, 3] = .15
capv0[4, 4] = .6
capv0[5, 5] = .6
capv0inv = solve(capv0)

{
# do posterior analysis
source("exercise_2_post.R")

# Posterior results based on Informative Prior
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
ystarm
ystarsd
ystarcapv
lmarglik
}

{
# Non-informative prior
v0 = 0
capv0inv = 0 * diag(k)

# do posterior analysis
source("exercise_2_post.R")

# Posterior results based on Non-informative Prior
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
ystarm
ystarsd
ystarcapv
}

# Likelihood tuning
marg_likelihood = vector("numeric")

b0 = c(0, 10, 5000, 10000, 10000)
s02 = 1 / 4.0e-8
capv0 = 2.4 * diag(k)
capv0[2, 2] = 6e-7
capv0[3, 3] = .15
capv0[4, 4] = .6
capv0[5, 5] = .6
capv0inv = solve(capv0)

alpha = c(0.0001, 0.1, 0.5, 1, 2, 10)
for(j in 1:length(alpha)) {
  v0 = alpha[j]
  source("exercise_2_post.R")
  marg_likelihood[j] = lmarglik
}
