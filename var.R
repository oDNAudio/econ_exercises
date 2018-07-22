bvar <- function(data,
                 lags,
                 constant = TRUE,
                 irf = TRUE, horizon = 20, nburn = 1000, nsave = 1000,
                 mn_prior = TRUE,
                 coeff_prior = FALSE,
                 no_coint = FALSE) {

  Y <- as.matrix(data)
  X <- lag_data(Y, lags)

  Y <- Y[(lags + 1):nrow(Y), ]
  X <- X[(lags + 1):nrow(X), ]


  if(constant) X <- cbind(X, 1)

  K <- ncol(X)
  M <- ncol(Y)
  N <- nrow(Y)

  # OLS
  B_ols <- solve(crossprod(X)) %*% crossprod(X, Y)
  SSE <- crossprod(Y- X %*% B_ols)
  sig_ols <- sig <- SSE / (N - K)

  ### Hyperpriors
  n_priors <- 3

  Y0 <- colMeans(Y[1:lags, ])

  par_min <- list("lambda" = 0.0001, "miu" = 0.0001, "theta" = 0.0001)
  par_max <- list("lambda" = 5, "miu" = 50, "theta" = 50)
  scale_c <- 0.1
  par_modes <- list("lambda" = 0.2, "miu" = 1, "theta" = 1)
  par_sd <- list("lambda" = 0.4, "miu" = 1, "theta" = 1)
  mn_alpha <- 2

  # Coefficients of the hyperpriors
  prior_coef <- list("lambda" = gamma_coef(par_modes$lambda, par_sd$lambda),
                     "miu" = gamma_coef(par_modes$miu, par_sd$miu),
                     "theta" = gamma_coef(par_modes$theta, par_sd$theta))

  # Inverse Hessian and Jacobian for drawing proposals
  H <- diag(n_priors) * 10
  exp_modes <- exp(c(par_modes$lambda, par_modes$miu, par_modes$theta))
  J <- exp_modes / (1 + exp_modes) ^ 2
  J[1] <- J[1] * (par_max$lambda - par_min$lambda)
  J[2] <- J[2] * (par_max$miu - par_min$miu)
  J[3] <- J[3] * (par_max$theta - par_min$theta)
  J <- diag(J)
  HH <- J %*% H %*% t(J)

  mn_mean <- matrix(0, K, M)
  mn_mean[2:(M + 1), ] <- diag(M)
  mn_sd <- apply(Y, 2, function(x) {
    sqrt(arima(x, order = c(lags, 0, 0))$sigma2)
  })
  mn_var <- 1e06

  post_mode <- c("lambda" = par_modes$lambda,
                 "miu" = par_modes$miu,
                 "theta" = par_modes$theta)

  in_bounds = FALSE
  while(!in_bounds) {
    par_draw <- MASS::mvrnorm(n = 1, mu = post_mode, Sigma = HH)
    if(!any(par_draw[1] < par_min[1] || par_draw[1] > par_max[1],
            par_draw[2] < par_min[2] || par_draw[2] > par_max[2],
            par_draw[3] < par_min[3] || par_draw[3] > par_max[3])) {
      in_bounds <- TRUE
    }
  }

  logML <- logML(Y, X, lags, par = par_draw, Y_row = N, Y_col = M,
                 mn_mean, mn_sd, mn_var, Y0, prior_coef)

  # # Add priors
  # if(mn_prior) {
  #   mn_prior <- minnesota_prior(Y, theta, gamma, delta, M, lags)
  #   Y <- rbind(m_prior$Y, Y)
  #   X <- rbind(m_prior$X, X)
  # }
  # if(coeff_prior) sum_coeff_prior()
  # if(prior3) prior3()

  # Get posteriors
  V_post <- solve(crossprod(X))
  A_post <- V_post %*% crossprod(X, Y)
  S_post <- crossprod(Y - X %*% A_post)
  S_post_inv <- solve(S_post)
  v_post <- nrow(Y)

  # Pre-calculate 1-step predictive density
  if(constant) {
    if(lags == 1) {
      X_pred <- c(Y[N, ], 1)
    } else {
      X_pred <- c(Y[N, ], X[N, 1:(M * (lags - 1))], 1)
    }
  } else {
    if(lags == 1) {
      X_pred <- Y[N, ]
    } else {
      X_pred <- c(Y[N, ], X[N, 1:(M * (lags - 1))])
    }
  }

  # Loop
  A_store <- array(0, c(nsave, K, M))
  Y_store <- matrix(0, nsave, M)
  IRF_store <- array(NA, c(nsave, M, M, horizon))

  for(i in 1:(nburn + nsave)) {
    # Posterior of A
    V_sim <- kronecker(sig, V_post)
    V_chol <- t(chol(V_sim))
    A_draw <- as.vector(A_post) + V_chol %*% rnorm(K * M)
    A_draw <- matrix(A_draw, K, M)
    # Posterior of Sigma
    sig_inv <- matrix(rWishart(1, v_post, S_post_inv), M, M)
    sig <- solve(sig_inv)

    if(i > nburn) {
      A_store[i - nburn, , ] <- A_draw

      # 1-step predicitve density
      Y_store[i - nburn, ] <- X_pred %*% A_draw + t(t(chol(sig)) %*% rnorm(M))

      ### Impulse responses

      # 1. Companion Matrix
      comp_size <- ifelse(constant, K - 1, K)
      B_comp <- matrix(0, comp_size, comp_size)
      B_comp[1:M, ] <- t(A_draw[1:comp_size, ])
      if(lags > 1) B_comp[(M + 1):comp_size, 1:(comp_size - M)] <- diag(M * (lags - 1))

      # 2. Compute IR
      shock <- t(chol(sig))

      IRF_draw <- array(0, c(M * lags, M * lags, horizon))
      IRF_draw[1:M, 1:M, 1] <- shock
      for(j in 2:horizon) {
        IRF_draw[, , j] <- IRF_draw[, , j - 1] %*% t(B_comp)
      }
      IRF_draw <- IRF_draw[1:M, 1:M, ]
      IRF_store[i - nburn, , ,] <- IRF_draw
    }
  }

  out <- list("A" = A_store, "Y" = Y_store, "IRF" = IRF_store)
}

lag_data <- function(x, lags) {

  x <- as.matrix(x)
  x_rows <- nrow(x)
  x_cols <- ncol(x)

  x_lag <- matrix(0, x_rows, lags * x_cols)

  for (i in 1:lags) {
    x_lag[(lags + 1):x_rows, (x_cols * (i - 1) + 1):(x_cols * i)] <-
      x[(lags + 1 - i):(x_rows - i), (1:x_cols)]
  }
  return(x_lag)
}

minnesota_prior <- function(x,
                            theta, gamma, delta,
                            lags) {

  sigma <- apply(x, 2, function(x) {
    sqrt(arima(x, order = c(lags, 0, 0))$sigma2)
  })
  M <- length(sigma)

  Y <- rbind(diag((as.numeric(sigma) * delta) / theta),
             matrix(0, M * (lags - 1), M),
             diag(as.numeric(sigma)),
             matrix(0, 1, M))

  X <- rbind(cbind(kronecker(diag(1:lags), diag(as.numeric(sigma))) / theta,
                   matrix(0, M * lags, 1)),
             matrix(0, M, M * lags + 1),
             cbind(matrix(0, 1, M * lags), gamma))

  return(list("Y" = Y, "X" = X))
}

gamma_coef <- function(mode, sd) {

  k <- (2 + mode^2 / sd^2 + sqrt((4 + mode^2 / sd^2) * mode^2 / sd^2)) / 2
  theta <- sqrt(sd ^ 2 / k)

  return(list("k" = k, "theta" = theta))
}
