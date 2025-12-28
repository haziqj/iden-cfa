library(lavaan)
library(tidyverse)

unit_to_stereo <- function(x) {
  # Project x \in R^p (a point on a sphere S^{p-1}) to y \in R^{p-1} via
  # stereographic projection from the north pole (0,...,0,1) onto the hyperplane
  # x_p = 0
  p <- length(x)
  out <- x[1:(p - 1)] / (1 - x[p])
  out
}

stereo_to_unit <- function(y) {
  # Inverse stereographic projection R^{p-1} -> S^{p-1} (unit sphere in R^{p})
  ysqnorm <- sum(y * y)
  out <- c(2 * y, ysqnorm - 1) / (1 + ysqnorm)
  out
}

source("R/partable.R")

geom_cfa <- function(
    model,
    data,
    method = c(""),
    start = NULL,
    control = list(),
    verbose = TRUE,
    add_priors = TRUE,
    ...
) {
  
  # Check arguments ------------------------------------------------------------
  method <- match.arg(method)
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$verbose <- verbose
  lavargs$do.fit <- FALSE
  
  # To prevent anchoring or fixing variance
  lavargs$auto.fix.first <- FALSE
  orig_std.lv <- lavargs$std.lv
  if (is.null(orig_std.lv)) orig_std.lv <- FALSE
  lavargs$std.lv <- FALSE
  
  # Initialise lavaan object ---------------------------------------------------
  fit0 <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options
  lavpartable    <- fit0@ParTable
  lavcache       <- fit0@Cache
  m              <- sum(lavpartable$free > 0)
  
  # Sample stats
  S <- fit0@SampleStats@cov
  n <- fit0@SampleStats@ntotal
  
  
  
  pt <- inlavaanify_partable(lavpartable, blavaan::dpriors(), lavdata, lavoptions)
  
  Sigma_from <- function(par) {
    # unpack
    z1 <- par[ 1:4]                 # spherical coords for f1 loadings (→ 5 entries)
    z2 <- par[ 5:8]                 # spherical coords for f2 loadings
    lv  <- par[ 9:10]               # log latent variances
    zr  <- par[11]                  # atanh latent correlation
    lt  <- par[12:21]               # log residual variances (10 items)
    
    L  <- Lambda_from(z1, z2)
    psi <- exp(lv)
    rho <- tanh(zr)
    
    Dhalf <- diag(sqrt(psi), 2, 2)
    Phi   <- matrix(c(1, rho, rho, 1), 2, 2)
    Psi   <- Dhalf %*% Phi %*% Dhalf
    
    Theta <- diag(exp(lt), 10, 10)
    L %*% Psi %*% t(L) + Theta
  }
  
  loglik_cov <- function(S, n, Sigma) {
    # -0.5 * n * (log|Σ| + tr(S Σ^{-1}))
    L <- try(chol(Sigma), silent = TRUE)
    if (inherits(L, "try-error")) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    # tr(S Σ^{-1}) via solves
    X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
    tr_term <- sum(diag(X))
    -0.5 * n * (logdet + tr_term)
  }
  
}

mod <- "
  f1 =~ x1 + x2 + x3
  f2 =~ x4 + x5 + x6
"

fit0 <- geom_cfa(mod, HolzingerSwineford1939)
















spherical_cfa <- function(S, n) {

  # Build Λ (p x 2) with simple-structure blocks (1:5 for f1, 6:10 for f2)
  Lambda_from <- function(z1, z2) {
    lam1 <- stereo_to_unit(z1)              # length 5, unit norm
    lam2 <- stereo_to_unit(z2)              # length 5, unit norm
    L <- matrix(0, nrow = 10, ncol = 2)
    L[1:5, 1]  <- lam1
    L[6:10, 2] <- lam2
    # fix orientation (optional): make first loading of each factor positive
    if (L[1,1] < 0)  L[1:5,1]  <- -L[1:5,1]
    if (L[6,2] < 0)  L[6:10,2] <- -L[6:10,2]
    L
  }
  
  # Build Σ = Λ Ψ Λ' + Θ from unconstrained params
  Sigma_from <- function(par) {
    # unpack
    z1 <- par[ 1:4]                 # spherical coords for f1 loadings (→ 5 entries)
    z2 <- par[ 5:8]                 # spherical coords for f2 loadings
    lv  <- par[ 9:10]               # log latent variances
    zr  <- par[11]                  # atanh latent correlation
    lt  <- par[12:21]               # log residual variances (10 items)
    
    L  <- Lambda_from(z1, z2)
    psi <- exp(lv)
    rho <- tanh(zr)
    
    Dhalf <- diag(sqrt(psi), 2, 2)
    Phi   <- matrix(c(1, rho, rho, 1), 2, 2)
    Psi   <- Dhalf %*% Phi %*% Dhalf
    
    Theta <- diag(exp(lt), 10, 10)
    L %*% Psi %*% t(L) + Theta
  }
  
  # Gaussian log-likelihood given sample covariance S and n
  loglik_cov <- function(S, n, Sigma) {
    # -0.5 * n * (log|Σ| + tr(S Σ^{-1}))
    L <- try(chol(Sigma), silent = TRUE)
    if (inherits(L, "try-error")) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    # tr(S Σ^{-1}) via solves
    X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
    tr_term <- sum(diag(X))
    -0.5 * n * (logdet + tr_term)
  }
  
  ## --- Objective -------------------------------------------------------------
  
  make_objective <- function(S, n) {
    force(S); force(n)
    function(par) {
      Sigma <- Sigma_from(par)
      # (Optional) weakly informative priors for stability:
      # normal(0,1) on z1,z2; N(0,1) on log-variances; N(0,1) on atanh(rho)
      lp <- sum(dnorm(par[1:8], 0, 1, log=TRUE)) +
        sum(dnorm(par[9:10], 0, 1, log=TRUE)) +
        dnorm(par[11], 0, 1, log=TRUE) +
        sum(dnorm(par[12:21], 0, 1, log=TRUE))
      # return NEGATIVE joint log-posterior (for minimisers)
      - (loglik_cov(S, n, Sigma) + lp)
    }
  }
  
  ## --- Initial values --------------------------------------------------------
  
  init_params <- function(S) {
    # crude starts: small angles (near north pole), unit latent vars, rho≈0,
    # residuals ~ log of diag minus small constant
    z1 <- rep(0, 4)
    z2 <- rep(0, 4)
    lv <- c(0, 0)                # log psi ~ 0 → psi ~ 1
    zr <- atanh(0.3)             # start rho at 0.3
    diagS <- pmax(1e-3, diag(S))
    lt <- log(pmax(1e-2, diagS * 0.3))
    c(z1, z2, lv, zr, lt)
  }
  
  ## --- Fit -------------------------------------------------------------------
  
  fit_spherical_cfa <- function(S, n) {
    fn <- make_objective(S, n)
    par0 <- init_params(S)
    opt <- optim(par0, fn, method = "BFGS", control = list(reltol = 1e-9, maxit = 5000))
    H   <- hessian(fn, opt$par)   # Hessian in unconstrained (spherical) coords
    list(par = opt$par, value = opt$value, conv = opt$convergence, hess = H)
  }
  
  ## --- Run on your data ------------------------------------------------------
  
  # Your S, n from the prompt:
  # S <- cov(dat); n <- nrow(dat)
  
  res <- fit_spherical_cfa(S, n)
  
  # Extract estimates
  par_hat <- res$par
  L_hat   <- Lambda_from(par_hat[1:4], par_hat[5:8])
  psi_hat <- exp(par_hat[9:10])
  rho_hat <- tanh(par_hat[11])
  Theta_hat <- diag(exp(par_hat[12:21]), 10, 10)
  
  # Quick checks
  c(norm_L1 = sqrt(sum(L_hat[1:5,1]^2)),
    norm_L2 = sqrt(sum(L_hat[6:10,2]^2)),
    psi_hat = psi_hat,
    rho_hat = rho_hat)
  
  # Condition number of on-manifold Hessian (scale-invariant)
  H  <- 0.5 * (res$hess + t(res$hess))
  eps <- sqrt(.Machine$double.eps)
  d  <- sqrt(pmax(eps, abs(diag(H))))
  S_H <- H / (d %o% d)
  ev  <- eigen(S_H, symmetric = TRUE, only.values = TRUE)$values
  ev  <- ev[ev > eps * max(ev)]
  kappa_H <- max(ev) / min(ev)
  kappa_H}
