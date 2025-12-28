library(tidyverse)
library(furrr)
library(pracma)   # for nullspace()
plan(multisession, workers = 14)

set.seed(123)
mod <- "
  f =~ y1 + y2 + y3 + y4 + y5
  y1 ~~ 0.1*y1  # more reliable
  y2 ~~ 0.2*y2
  y3 ~~ 0.3*y3
  y4 ~~ 0.4*y4
  y5 ~~ 1*y5   # less reliable
"
n <- 25000
dat <- lavaan::simulateData(mod, sample.nobs = n)
S <- cov(dat)
p <- ncol(S)

# mean-zero orthonormal basis H (p x (p-1)) with 1' H = 0 and H'H = I
mean_zero_basis <- function(p) {
  one <- rep(1, p) / sqrt(p)
  Q <- qr.Q(qr(cbind(one, diag(p))))[, -1, drop = FALSE]
  Q
}
H <- mean_zero_basis(p)  # fixed once

# stable zero-mean MVN log-lik from sufficient stats S
dmvnorm0 <- function(Sigma, S, n) {
  L <- try(chol(Sigma), silent = TRUE)
  if (inherits(L, "try-error")) return(-Inf)
  logdet <- 2 * sum(log(diag(L)))
  X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
  tr_term <- sum(diag(X))
  -0.5 * n * (logdet + tr_term)
}

# Planar (effects-coded) profile log-likelihood over (lambda_i, log psi)
# Lambda = 1 + H b, with mean(lambda)=1 enforced intrinsically.
# Fix lambda_i on the grid by solving H[i,] b = lambda_i - 1 -> b = b0 + N u, N spans null(H[i,])
loglik_planar <- function(lambda_i, logpsi, i = 1) {
  
  psi <- exp(logpsi)
  if (!is.finite(psi) || psi <= 0) return(-Inf)
  
  # Row of H corresponding to item i (length p-1)
  h_i <- H[i, , drop = FALSE]               # 1 x (p-1)
  rhs <- lambda_i - 1
  
  # Minimal-norm solution b0 to h_i b = rhs
  denom <- as.numeric(h_i %*% t(h_i))       # ||h_i||^2
  if (denom < .Machine$double.eps) return(-Inf)
  b0 <- as.numeric(t(h_i) * (rhs / denom))  # (p-1) vector
  
  # Nullspace N of the row constraint (dimension p-2)
  N <- nullspace(h_i)                       # (p-1) x (p-2)
  d_free <- if (is.null(dim(N))) 0L else ncol(N)
  
  fn <- function(par) {
    log_theta <- par[seq_len(p)]
    theta <- exp(log_theta)
    if (any(!is.finite(theta)) || any(theta <= 0)) return(Inf)
    
    if (d_free > 0) {
      u <- par[(p + 1):(p + d_free)]
      b <- b0 + as.vector(N %*% u)
    } else {
      b <- b0
    }
    
    lambda <- as.vector(1 + H %*% b)       # effects-coded loadings (mean = 1)
    Lambda <- matrix(lambda, ncol = 1)
    Sigma  <- psi * tcrossprod(Lambda) + diag(theta)
    # return NEGATIVE log posterior (no priors here; add if you like)
    -dmvnorm0(Sigma, S, n)
  }
  
  # starts: log theta from diag(S); u = 0
  par0 <- c(
    log(pmax(diag(S) * 0.5, 1e-3)),
    if (d_free > 0) rep(0, d_free)
  )
  
  opt <- nlminb(start = par0, objective = fn,
                control = list(rel.tol = 1e-8, x.tol = 1e-8,
                               eval.max = 1000, iter.max = 1000))
  if (is.finite(opt$objective)) -opt$objective else -Inf
}

# Evaluate true profile (lambda_i, log psi) with i = high vs low reliability
plot_df <-
  expand_grid(
    lambda = seq(0.01, 4, length.out = 250),
    logpsi = seq(-6,   6, length.out = 250),
    i = c(1, 5)
  ) |>
  mutate(
    ll = future_pmap_dbl(
      .l = list(lambda = lambda, logpsi = logpsi, i = i),
      .f = loglik_planar,
      .progress = TRUE
    ),
    i = factor(i, labels = c("High reliability", "Low reliability"))
  )

p_profll_planar <-
  ggplot(plot_df, aes(lambda, logpsi)) +
  geom_raster(aes(fill = ll)) +
  geom_contour(aes(z = ll), col = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 1,     linetype = "dashed", color = "black") +
  scale_fill_viridis_c(option = "turbo", labels = scales::label_comma()) +
  labs(
    x = expression(lambda[i]),
    y = expression(log(psi)),
    fill = "Log-likelihood    "
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width  = unit(2, "cm"),
    legend.key.height = unit(0.2, "cm")
  ) +
  facet_grid(. ~ i)

p_profll_planar
save(p_profll_planar, file = "profll_planar.RData")