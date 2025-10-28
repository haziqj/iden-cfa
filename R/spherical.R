library(tidyverse)
library(furrr)
plan(multisession, workers = 14)

set.seed(123)
mod <- "
  f =~ y1 + y2 + y3 + y4 + y5
  y1 ~~ 0.1*y1  # more reliable
  y2 ~~ 0.2*y2
  y3 ~~ 0.3*y3
  y4 ~~ 0.4*y4
  y5 ~~ 1*y5  # less reliable
"
n <- 25000
dat <- lavaan::simulateData(mod, sample.nobs = n)
S <- cov(dat)
p <- ncol(S)

# build unit direction tildeLambda with fixed coordinate i set to tl_i in (-1,1)
# a: unconstrained R^{p-1} params for the remaining components (we normalise them)
make_unit_direction <- function(a, tl_i, i, p) {
  # handle zero vector 'a' (use equal weights)
  if (all(!is.finite(a)) || sum(a^2) == 0) {
    a <- rep(1, p - 1)
  }
  a <- a / sqrt(sum(a ^ 2))                 # unit vector in R^{p-1}
  scale_rest <- sqrt(max(1 - tl_i ^ 2, 0))  # enforce ||tildeLambda|| = 1
  tl <- numeric(p)
  tl[i] <- tl_i
  tl[-i] <- scale_rest * a                  # fill the remaining coords
  tl
}

# stable zero-mean MVN log-lik from sufficient stats S (your dmvnorm0)
dmvnorm0 <- function(Sigma, S, n) {
  L <- try(chol(Sigma), silent = TRUE)
  if (inherits(L, "try-error")) return(-Inf)
  logdet <- 2 * sum(log(diag(L)))
  X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
  tr_term <- sum(diag(X))
  -0.5 * n * (logdet + tr_term)
}

loglik_spherical <- function(lambda_i, logpsi, i = 1) {
  
  r   <- exp(0.5 * logpsi)             # scale = sqrt(psi)
  if (!is.finite(r) || r <= 0) return(-Inf)
  
  tl_i <- lambda_i / r                  # implied standardized loading coord
  if (!is.finite(tl_i) || abs(tl_i) >= 1) return(-Inf)  # infeasible region on sphere
  
  # objective over nuisance: log(theta_1..theta_p) and a (p-1 dims for direction)
  # parameter vector: c(log_theta[1:p], a[1:(p-1)])
  fn <- function(par) {
    log_theta <- par[1:p]
    a         <- par[(p+1):(p + (p-1))]
    theta     <- exp(log_theta)
    
    tl <- make_unit_direction(a, tl_i, i, p)     # tildeLambda (unit)
    Lambda <- r * tl                             # Lambda = r * tildeLambda
    Sigma  <- tcrossprod(Lambda) + diag(theta) * 1.0  # r^2 already in tcrossprod
    
    # NOTE: tcrossprod(Lambda) = (r*tl)(r*tl)' = r^2 * tl tl'
    dmvnorm0(Sigma, S, n)
  }
  
  # starting values (robust & cheap)
  par0 <- c(
    log(pmax(diag(S) * 0.5, 1e-3)),     # log theta start
    rep(1, p - 1)                       # a start (direction filler)
  )
  
  opt <- nlminb(start = par0, objective = function(z) -fn(z),
                control = list(rel.tol = 1e-8, x.tol = 1e-8, eval.max = 1000, iter.max = 1000))
  if (is.finite(opt$objective)) return(-opt$objective) else return(-Inf)
}

plot_df <-
  expand_grid( 
    lambda = seq(0.01, 4, length.out = 250),
    logpsi = seq(-6, 6, length.out = 250),
    i = c(1, 5)
  ) |>
  mutate(
    ll = future_pmap_dbl(
      .l = list(lambda = lambda, logpsi = logpsi, i = i),
      .f = loglik_spherical,
      .progress = TRUE),
    i = factor(i, labels = c("High reliability", "Low reliability"))
  ) 

p_profll_spherical <-
  ggplot(plot_df, aes(lambda, logpsi)) +
  geom_raster(aes(fill = ll)) +
  geom_contour(aes(z = ll), col = "white") +
  geom_hline(yintercept = log(5), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_viridis_c(option = "turbo",
    labels = scales::label_comma()
  ) +
  # scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(expand = c(0, 0)) +
  labs(
    # title = expression("Profile log-likelihood surface over ("*lambda[i]*", "*psi*")"),
    x = expression(lambda[i]),
    y = expression(log(psi)),
    fill = "Log-likelihood    "
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.2, "cm")
  ) +  
  facet_grid(. ~ i); p_profll_spherical

save(p_profll_spherical, file = "profll_spherical.RData")
