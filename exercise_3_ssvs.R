#' @title Stochastic Search Variable Selection
#' @author Nikolas Kuschnig
#' @description Uses the Gibbs sampler to perform SSVS
#'
#' @param y The endogenous variable, must be convertible to a matrix.
#' @param X The explanatory variables, must be convertible to a matrix.
#' @param tau0 Number by which to scale an "unimportant" variable. Will use least squares estimates if not supplied.
#' @param tau1 Number by which to scale an "important" variable. Defaults to 100 * tau0.
#' @param save Iterations to be stored for calculating means.
#' @param burn Iterations to be discarded before calculating means.
#' @param s_prior Prior accuracy, i.e. weight assigned to the prior
#' @param S_prior 1 / sigma^2
#' @param standardise A boolean determining whether to center and scale X & y.
#'
#' @return Returns a list containing the means of posterior: inclusion probability, mean and standard deviation.
#' @export

ssvs = function(y, 
                X, 
                tau0 = NULL, 
                tau1 = tau0 * 100, 
                save = 4000, 
                burn = 1000, 
                s_prior = 0.01,
                S_prior = 0.01,
                standardise = TRUE) {
  y = matrix(y)
  X = as.matrix(X)

  N = nrow(y)
  K = ncol(X)
  
  # Standardise (i.e. center and scale) the data if desired
  if(standardise) {
    y = scale(y)
    X = scale(X)
  }
  
  # If no tau0 was supplied set up for OLS estimates, otherwise vectorise
  if(is.null(tau0)) {
    c0 = 0.1
    c1 = ifelse(length(tau1) == 0, 10, tau1)
  } else {
    tau0 = rep(tau0, K)
    tau1 = rep(tau1, K)
  }
  
  # get OLS stuff
  OLS = solve(crossprod(X)) %*% crossprod(X, y)
  SSE = as.numeric(crossprod(y - X %*% OLS))
  sigma_draw = as.numeric(SSE / (N - K))
  V_beta_draw = sigma_draw * solve(crossprod(X))
  
  # use OLS estimates for tau if no tau0 was supplied
  if(is.null(tau0)) {
    tau0 <- c0 * sqrt(diag(V_beta_draw))
    tau1 <- c1 * sqrt(diag(V_beta_draw))
  }
  
  # Set up gamma, and storage matrices
  gamma = matrix(1, K, 1)
  V_prior = diag(as.numeric(gamma * tau1 + (1 - gamma) * tau0))

  alpha_store = matrix(NA, save, K)
  sigma_store = matrix(NA, save, 1)
  V_beta_store = matrix(NA, save, K)
  gamma_store = matrix(NA, save, K)
  
  # Do the loop
  for(i in 1:(save + burn)) {
    # Draw alpha
    V_post = solve(crossprod(X) / sigma_draw + diag(1 / diag(V_prior)))
    alpha_post = V_post %*% (crossprod(X, y) / sigma_draw)
    alpha_draw = alpha_post + t(chol(V_post)) %*% rnorm(K)
    
    # Determine inclusion based on alpha
    for(j in 1:K) {
      p0 = dnorm(alpha_draw[[j]], 0, sqrt(tau0[j]))
      p1 = dnorm(alpha_draw[[j]], 0, sqrt(tau1[j]))
      p11 = p1 / (p0 + p1)
      
      gamma[[j]] = ifelse(p11 > runif(1), 1, 0)
    }
    
    # Construct prior VC matrix
    V_prior = diag(as.numeric(gamma * tau1 + (1 - gamma) * tau0))
    
    # Draw sigma^2
    S_post = S_prior + crossprod(y - X %*% alpha_draw) / 2
    s_post = S_prior + N / 2
    sigma_draw = 1 / rgamma(1, s_post, S_post)
    V_beta_draw = diag(sigma_draw * solve(crossprod(X)))
    
    # Ignore the first n(=burn) iterations, store results afterwards
    if(i > burn) {
      alpha_store[i - burn, ] = alpha_draw
      sigma_store[i - burn, ] = sigma_draw
      V_beta_store[i - burn, ] = V_beta_draw
      gamma_store[i - burn, ] = gamma
    }
  }
  
  # Get means for the output
  pip_mean = apply(gamma_store, 2, mean)
  alpha_mean = apply(alpha_store, 2, mean)
  sigma_mean = apply(sigma_store, 2, mean)
  V_beta_mean = apply(V_beta_store, 2, mean)
  
  names(pip_mean) = names(alpha_mean) = names(V_beta_mean) = colnames(X)
  
  # Store the parameters used for the output
  meta_data = list("tau0" = tau0, "tau1" = tau1,
                   c("save" = save, "burn" = burn,
                   "s_prior" = s_prior, "S_prior" = S_prior,
                   "standardise" = standardise))
  
  out = list(pip_mean, alpha_mean, sigma_mean, V_beta_mean, meta_data)
  
  names(out) = c("pip", "post_mean", "post_std", "post_var_beta", "meta")

  return(out)
}

