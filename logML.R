logML <- function(Y, X, lags, par, Y_row, Y_col, min = NULL, max = NULL, 
                  mn_mean, mn_sd, mn_var, Y0, prior_coef) {
    
    if(!is.null(min) && !is.null(max)) {
        if(any(par$lambda < min$lambda || par$lambda < max$lambda, 
            par$psi < min$psi || par$psi < max$psi, 
            par$theta < min$theta || par$theta < max$theta, 
            par$miu < min$miu || par$miu < max$miu, 
            par$alpha < min$alpha || par$alpha < max$alpha)) {

            return(list(logML = -10e15, beta_draw = NULL, sigma_draw = NULL))
            #stop("Parameters outside bounds - Quitting.")
        }
    }
    
    omega <- vector("numeric", 1 + Y_col * lags)
    omega[1] <- mn_var
    for(i in 1:lags) {
        omega[seq(2 + Y_col * (i - 1), 1 + i * Y_col)] <- (par$lambda^2 / i^par$alpha) / par$psi
    }
    psi <- diag(par$psi)

    # single-unit-root & sum-of-coefficients
    sur <- Y0 / par$theta
    soc <- diag(Y0)
    Y <- rbind(Y, 
               sur, 
               soc / par$miu)
    X <- rbind(X, 
               c(1 / par$theta, rep(sur, lags)), 
               cbind(rep(0, Y_row), rep(soc, lags) / par$miu))
    
    Dummy_row <- nrow(Y) - Y_row
    Y_row <- nrow(Y)
    
    omega_inv <- diag(1 / omega)
    beta_hat <- solve(crossprod(X) + omega_inv) %*% (crossprod(X, Y) + omega_inv * mn_mean)
    sse <- crossprod(Y - X %*% beta_hat)

    psi_inv <- diag(1 / sqrt(psi))
    aaa <- diag(sqrt(omega)) %*% crossprod(X) %*% diag(sqrt(omega))
    bbb <- psi_inv %*% (sse + t(beta_hat - mn_mean) %*% omega_inv %*% (beta_hat - mn_mean)) %*% psi_inv

    aaa_eig <- Re(eigen(aaa, only.values = TRUE))
    aaa_eig[aaa_eig < 1e-12] <- 0
    aaa_eig <- aaa_eig + 1
    bbb_eig <- Re(eigen(bbb, only.values = TRUE))
    bbb_eig[bbb_eig < 1e-12] <- 0
    bbb_eig <- bbb_eig + 1

    logML <- (-Y_col * Y_row * log(pi) / 2) + 
        sum(lgamma(((Y_row + Y_col + 2) - 0:(Y_col - 1)) / 2)) - 
        lgamma(((Y-col + 2) - 0:(Y_col -1)) / 2) - 
        (Y_row * sum(log(psi)) / 2) - 
        (Y_col * sum(log(aaa_eig)) / 2) + 
        ((Y_row + Y_col + 2) * sum(log(bbb_eig)) / 2)

    # here we would normalise

    logML <- logML + 
        lgamma_pdf(par$lambda, prior_coef$lambda_k, prior_coef$lambda_theta) +
        lgamma_pdf(par$theta, prior_coef$theta_k, prior_coef$theta_theta) +
        lgamma_pdf(par$miu, prior_coef$miu_k, prior_coef$miu_theta)

    S <- psi + sse + t(beta_hat - mn_mean) %*% omega_inv %*% (beta_hat - mn_mean)
    S_eig <- eigen(S)
    S_inv <- S_eig$vectors %*% diag(1 / abs(diag(S_eig$values))) %*% t(S_eig$vectors)
    eta <- MASS:mvrnorm(n = (Y_row + Y_col + 2), mu = rep(0, Y_col), Sigma = S_inv)
    sigma_draw <- solve(crossprod(eta)) %*% diag(Y_col)
    sigma_chol <- chol(sigma_draw)
    beta_draw <- beta_hat + 
        mvrnorm(n = Y_col, 
            mu = rep(0, (1 + Y_col * lags)), 
            Sigma = crossprod(X) + omega_inv) %*% sigma_chol

    return(list(logML = logML, beta_draw = beta_draw, sigma_draw = sigma_draw))
}

lgamma_pdf <- function(x, k, theta) {
    r <- (k - 1) * log(x) - x / theta - k * log(theta) - lgamma(k)
}