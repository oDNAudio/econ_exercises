library(bvarsv)
source("aux_functions.R")

data(usmacro)


# data and configuration setup
p <- 4
cons <- TRUE
nhor <- 13 # number of predictions

Y <- as.matrix(usmacro)
X <- mlag(Y,p)
Y <- Y[(p+1):nrow(Y),]
X <- X[(p+1):nrow(X),]
if (cons){
  X <- cbind(1,X)
}

K <- ncol(X)
M <- ncol(Y)


# OLS quantities
B.OLS <- solve(crossprod(X))%*%crossprod(X,Y)
yfit <- X%*%B.OLS

SSE <- crossprod(Y-X%*%B.OLS)
T <- nrow(Y)
SIGMA <- SIGMA.OLS <- SSE/(T-K)

# number of priors; all included, thus equal to 3
n.prior <- 3

MIN <- list()
MIN$lambda <- 0.0001
MIN$miu <- 0.0001
MIN$theta <- 0.0001

MAX <- list()
MAX$lambda <- 5
MAX$miu <- 50
MAX$theta <- 50

scale.c <- 0.1

#  hyperpriors modes
mode.lambda <- .2
mode.miu <- 1
mode.theta <- 1
#  hyperpriors std
sd.lambda <- .4      
sd.miu <- 1
sd.theta <- 1
# scale and shape of the IG on psi/(d-n-1)    
#scalePSI <- 0.02^2  
  
 
  # coefficients of hyperpriors
priorcoef <- list()
priorcoef$lambda <- GammaCoef(mode.lambda, sd.lambda)
priorcoef$miu <- GammaCoef(mode.miu,sd.miu)
priorcoef$theta <- GammaCoef(mode.theta,sd.theta)

gamma.prior <- 10^6 #prior on the intercept of the Minnesota prior

# Guess for inverse Hessian
H <- diag(n.prior)*10

# Jacobian used as vcov matrix for drawing the proposals
init.guesses <- c(mode.lambda, mode.miu, mode.theta)
JJ <- exp(init.guesses)/(1+exp(init.guesses))^2
JJ[1] <- (MAX$lambda-MIN$lambda)*JJ[1]
JJ[2] <- (MAX$miu-MIN$miu)*JJ[2]
JJ[3] <- (MAX$theta-MIN$theta)*JJ[3]
JJ <- diag(as.vector(JJ))
HH <- JJ%*%H%*%t(JJ)



## Minnesota prior values
mn.mean <- matrix(0, K, M)
mn.mean[2:(M+1),] <- diag(M)

mn.sd <- matrix(NA,M,1)
for (ii in 1:M){
  tmpar <- arima(Y[,ii],order=c(p,0,0))
  mn.sd[ii,1] <- sqrt(tmpar$sigma2)
}

mn.var <- 1e6
 
### starting values for Metropolis-Hastings
# 
# lambda.draw <- rnorm(1,init.guesses[1],diag(HH)[1])
# mu.draw <- rnorm(1,init.guesses[2],diag(HH)[2])
# delta.draw <- rnorm(1,init.guesses[3],diag(HH)[3]) 

draw <- list()
draw$lambda <- mode.min
draw$miu <- mode.soc
draw$theta <- mode.ido
draw$psi <- mn.sd
draw$alpha <- 2
prop <- draw


Y0 <- colMeans(Y[1:p,])



logML.old <- logML(Y=Y, X=X, lags=p, par=draw, Y_row=T, Y_col=M, mn_mean=mn.mean, mn_sd=mn.sd, mn_var=mn.var, Y0=Y0, prior_coef=priorcoef, min = MIN, max=MAX)



# get starting value for ML via function





nsave <- 5000
nburn <- 10000
ntot <- nsave+nburn
beta.store <- array(0,c(nsave,K,M))
Sigma.store <- array(0,c(nsave,K,(K-1)))
hypara.store <- matrix(0,nsave,3)
logML.store <- matrix(NA,nsave,1)

for(irep in 1:ntot){
  
  
  
  
  
  prop$lambda <- rnorm(1,draw$lambda,diag(HH)[1]*scale.c^2)
  prop$miu <- rnorm(1,draw$miu,diag(HH)[2]*scale.c^2)
  prop$theta <- rnorm(1,draw$theta,diag(HH)[3]*scale.c^2)
  
  ## get values for new ML via function
  temp <- logML(inputs)
  logML.prop <- temp$logML
  
  
  
  
  if((logML.prop-logML.old)>log(runif(1))){
    
    logML.old <- logML.prop
    draw <- prop
    
  }else{
    
    ## draw beta and Sigma with old values for hyperparameters
    # Funktion von Niko für neue Outputs
    
  }
  
  
  
  if(irep>nburn){
    beta.store[(irep-nburn),,] <- temp$beta_draw
    Sigma.store[(irep-nburn),,] <- temp$sigma_draw
    hypara.store[(irep-nburn),] <- c(draw$lamba, draw$miu, draw$theta)
    logML.store[(irep-nburn)] <- logML.old
  }
  
}