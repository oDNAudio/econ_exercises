ssvs = function(Y,
                    X,
                    tau0 = 0.01,
                    tau1 = tau0 * 100,
                    iterations = 1000,
                    burn = 1000,
                    s_prior = 0.01,
                    S_prior = 0.01,
                    standardise = FALSE) {
  Y = matrix(Y)
  X = as.matrix(cbind(rnorm(nrow(X), 0, 10), rnorm(nrow(X), 0, 1), X))

  if(standardise) {
    Y = scale(Y)
    X = scale(X)
  }
  
  N = nrow(Y)
  K = ncol(X)
  
  OLS = solve(crossprod(X)) %*% crossprod(X, Y)
  SSE = as.numeric(crossprod(Y - X %*% OLS))
  sigma_draw = as.numeric(SSE / (N - K))
  
  gamma = matrix(1, K, 1)

  alpha_store = matrix(NA, iterations, K)
  sigma_store = matrix(NA, iterations, 1)
  gamma_store = matrix(NA, iterations, K)
  
  for(i in 1:(iterations + burn)) {
    V_prior = diag(as.numeric(gamma * tau1 + (1 - gamma) * tau0))

    V_post = solve(crossprod(X) / sigma_draw + diag(1 / diag(V_prior)))
    A_post = V_post %*% (crossprod(X, Y) / sigma_draw)
    A_draw = A_post + t(chol(V_post)) %*% rnorm(K)
    
    for(j in 1:K) {
      p0 = dnorm(A_draw[[j]], 0, sqrt(tau0))
      p1 = dnorm(A_draw[[j]], 0, sqrt(tau1))
      p11 = p1 / (p0 + p1)
      
      gamma[[j]] = ifelse(p11 > runif(1), 1, 0)
    }
    
    S_post = S_prior + crossprod(Y - X %*% A_draw) / 2
    s_post = S_prior + N / 2
    
    sigma_draw = 1 / rgamma(1, s_post, S_post)
    
    if(i > burn) {
      alpha_store[i - burn, ] = A_draw
      sigma_store[i - burn, ] = sigma_draw
      gamma_store[i - burn, ] = gamma
    }
  }
  
  PIP_mean = apply(gamma_store, 2, mean)
  alpha_mean = apply(alpha_store, 2, mean)
  sigma_mean = apply(sigma_store, 2, mean)
  
  out = list(PIP_mean, alpha_mean, sigma_mean)
  #out = list(gamma_store, alpha_store, sigma_store)
  
  return(out)
}
