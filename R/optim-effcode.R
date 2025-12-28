planar_cfa <- function(S, n) {
  
  ## --- Mean-zero orthonormal basis (length p -> p x (p-1)) --------------------
  mean_zero_basis <- function(p) {
    one <- rep(1, p) / sqrt(p)
    Q <- qr.Q(qr(cbind(one, diag(p))))[, -1, drop = FALSE]  # columns ⟂ to 'one'
    Q  # p x (p-1), Q'Q = I, and t(one) %*% Q = 0
  }
  
  # Fixed basis for 5 items per factor
  H5 <- mean_zero_basis(5)  # 5 x 4
  
  ## --- Λ from intrinsic effect-coding coordinates -----------------------------
  # For each factor j, free params b_j ∈ R^{4}; loadings λ_j = 1 + H5 b_j  (mean(λ_j)=1)
  Lambda_from_ec <- function(b1, b2) {
    lam1 <- as.vector(1 + H5 %*% b1)  # length 5
    lam2 <- as.vector(1 + H5 %*% b2)
    L <- matrix(0, nrow = 10, ncol = 2)
    L[1:5, 1]  <- lam1
    L[6:10, 2] <- lam2
    L
  }
  
  ## --- Σ(θ) -------------------------------------------------------------------
  Sigma_from_ec <- function(par) {
    b1 <- par[ 1:4]
    b2 <- par[ 5:8]
    lv <- par[ 9:10]    # log latent variances
    zr <- par[11]       # Fisher-z latent correlation
    lt <- par[12:21]    # log residual variances (10 items)
    
    L   <- Lambda_from_ec(b1, b2)
    psi <- exp(lv)
    rho <- tanh(zr)
    
    Dhalf <- diag(sqrt(psi), 2, 2)
    Phi   <- matrix(c(1, rho, rho, 1), 2, 2)
    Psi   <- Dhalf %*% Phi %*% Dhalf
    
    Theta <- diag(exp(lt), 10, 10)
    L %*% Psi %*% t(L) + Theta
  }
  
  ## --- Gaussian log-likelihood on covariances --------------------------------
  loglik_cov <- function(S, n, Sigma) {
    L <- try(chol(Sigma), silent = TRUE)
    if (inherits(L, "try-error")) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
    tr_term <- sum(diag(X))
    -0.5 * n * (logdet + tr_term)
  }
  
  ## --- Objective: negative joint log-posterior (weakly informative priors) ----
  make_objective_ec <- function(S, n) {
    force(S); force(n)
    function(par) {
      Sigma <- Sigma_from_ec(par)
      lp <- sum(dnorm(par[1:8], 0, 1, log = TRUE)) +     # b1,b2
        sum(dnorm(par[9:10], 0, 1, log = TRUE)) +    # log psi
        dnorm(par[11], 0, 1, log = TRUE) +           # Fisher-z rho
        sum(dnorm(par[12:21], 0, 1, log = TRUE))     # log theta
      -(loglik_cov(S, n, Sigma) + lp)
    }
  }
  
  ## --- Inits ------------------------------------------------------------------
  init_params_ec <- function(S) {
    b1 <- rep(0, 4); b2 <- rep(0, 4)   # start near mean(λ)=1, no contrasts
    lv <- c(0, 0)                      # psi ≈ 1
    zr <- atanh(0.3)                   # rho ≈ 0.3
    diagS <- pmax(1e-3, diag(S))
    lt <- log(pmax(1e-2, diagS * 0.3)) # modest residuals
    c(b1, b2, lv, zr, lt)
  }
  
  ## --- Fit in intrinsic effect-coding coordinates -----------------------------
  fit_effectcoded_cfa <- function(S, n) {
    fn  <- make_objective_ec(S, n)
    p0  <- init_params_ec(S)
    opt <- optim(p0, fn, method = "BFGS",
                 control = list(reltol = 1e-9, maxit = 5000))
    H   <- hessian(fn, opt$par)
    list(par = opt$par, value = opt$value, conv = opt$convergence, hess = H)
  }
  
  ## === Run on your data =======================================================
  # S <- cov(dat); n <- nrow(dat)   # you already have these
  res_ec <- fit_effectcoded_cfa(S, n)
  
  # Extract estimates
  ph <- res_ec$par
  L_hat   <- Lambda_from_ec(ph[1:4], ph[5:8])
  psi_hat <- exp(ph[9:10])
  rho_hat <- tanh(ph[11])
  Theta_hat <- diag(exp(ph[12:21]), 10, 10)
  
  # Quick sanity
  c(mean_lambda_f1 = mean(L_hat[1:5,1]),
    mean_lambda_f2 = mean(L_hat[6:10,2]),
    psi_hat = psi_hat,
    rho_hat = rho_hat)
  
  # Condition number of intrinsic Hessian (scale-invariant)
  H  <- 0.5 * (res_ec$hess + t(res_ec$hess))
  eps <- sqrt(.Machine$double.eps)
  d  <- sqrt(pmax(eps, abs(diag(H))))
  S_H <- H / (d %o% d)
  ev  <- eigen(S_H, symmetric = TRUE, only.values = TRUE)$values
  ev  <- ev[ev > eps * max(ev)]
  kappa_H_ec <- max(ev) / min(ev)
  kappa_H_ec
}



