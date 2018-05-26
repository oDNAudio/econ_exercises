#' 
#' @title Stochastic Search Variable Selection
#' @author Nikolas Kuschnig
#' @description Uses the Gibbs sampler to perform SSVS
#'
#' @param Y The endogenous variable, must be convertible to a matrix.
#' @param X The explanatory variables, must be convertible to a matrix.
#' @param tau0 Number by which to scale an "unimportant" variable. Will use least squares estimate if not supplied.
#' @param tau1 Number by which to scale an "important" variable. Defaults to 100 * tau0.
#' @param save Iterations to be stored for calculating means.
#' @param burn Iterations to be discarded before calculating means.
#' @param s_prior Prior accuracy, i.e. weight assigned to the prior
#' @param S_prior 1 / sigma^2
#' @param standardise A boolean determining whether to center and scale X & Y.
#'
#' @return Returns a list containing the means of posterior: inclusion probability, mean and standard deviation.
#' @export

ssvs = function(Y, 
                X, 
                tau0 = NULL,
                tau1 = tau0 * 100, 
                save = 9000, 
                burn = 1000, 
                s_prior = 0.01,
                S_prior = 0.01,
                standardise = TRUE) {
  Y = matrix(Y)
  X = as.matrix(cbind(rnorm(nrow(X), 0, 10), rnorm(nrow(X), 0, 1), X))

  N = nrow(Y)
  K = ncol(X)
  
  if(standardise) {
    Y = scale(Y)
    X = scale(X)
  }
  
  if(is.null(tau0)) {
    c0 = 0.1
    c1 = 10
  } else {
    tau0 = rep(tau0, K)
    tau1 = rep(tau1, K)
  }
  
  OLS = solve(crossprod(X)) %*% crossprod(X, Y)
  SSE = as.numeric(crossprod(Y - X %*% OLS))
  sigma_draw = as.numeric(SSE / (N - K))
  V_beta = sigma_draw * solve(crossprod(X))
  
  if(is.null(tau0)) {
    tau0 <- c0 * sqrt(diag(V_beta))
    tau1 <- c1 * sqrt(diag(V_beta))
  }

  gamma = matrix(1, K, 1)
  V_prior = diag(as.numeric(gamma * tau1 + (1 - gamma) * tau0))

  alpha_store = matrix(NA, save, K)
  sigma_store = matrix(NA, save, 1)
  gamma_store = matrix(NA, save, K)
  
  for(i in 1:(save + burn)) {
    V_post = solve(crossprod(X) / sigma_draw + diag(1 / diag(V_prior)))
    alpha_post = V_post %*% (crossprod(X, Y) / sigma_draw)
    alpha_draw = alpha_post + t(chol(V_post)) %*% rnorm(K)
    
    for(j in 1:K) {
      p0 = dnorm(alpha_draw[[j]], 0, sqrt(tau0[j]))
      p1 = dnorm(alpha_draw[[j]], 0, sqrt(tau1[j]))
      p11 = p1 / (p0 + p1)
      
      gamma[[j]] = ifelse(p11 > runif(1), 1, 0)
    }
    
    V_prior = diag(as.numeric(gamma * tau1 + (1 - gamma) * tau0))
    
    S_post = S_prior + crossprod(Y - X %*% alpha_draw) / 2
    s_post = S_prior + N / 2
    
    sigma_draw = 1 / rgamma(1, s_post, S_post)
    
    if(i > burn) {
      alpha_store[i - burn, ] = alpha_draw
      sigma_store[i - burn, ] = sigma_draw
      gamma_store[i - burn, ] = gamma
    }
  }
  
  pip_mean = apply(gamma_store, 2, mean)
  alpha_mean = apply(alpha_store, 2, mean)
  sigma_mean = apply(sigma_store, 2, mean)
  
  out = list(pip_mean, alpha_mean, sigma_mean)
  
  names(out) = c("Posterior Inclusion Probability", "Posterior Mean", "Posterior Standard Deviation")

  return(out)
}

