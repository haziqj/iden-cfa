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

dmvnorm0 <- function(Sigma) {
  L <- try(chol(Sigma), silent = TRUE)
  if (inherits(L, "try-error")) return(-Inf)
  logdet <- 2 * sum(log(diag(L)))
  # Trace term via solves: tr(S Sigma^{-1}) = sum((L^{-1} S L^{-T}))
  X <- backsolve(L, t(backsolve(L, S, transpose = TRUE)), transpose = TRUE)
  tr_term <- sum(diag(X))
  -0.5 * n * (logdet + tr_term)
}

loglik <- function(lambda, logpsi, i) {
  
  psi <- exp(logpsi)
  
  # Profile the likelihood for loading i
  fn <- function(par) {
    theta <- exp(par[1:5])
    lambdavec <- numeric(5)
    lambdavec[-i] <- par[6:9]
    lambdavec[i] <- lambda
    Lambda <- matrix(lambdavec, ncol = 1)
    Sigmay <- psi * (Lambda %*% t(Lambda)) + diag(theta)
    -1 * dmvnorm0(Sigmay)
  }
  opt <- nlminb(c(log(c(0.1, 0.1, 0.3, 0.4, 0.5)), rep(1, 4)), fn)
  -1 * opt$objective
}

plot_df <-
  expand_grid( 
    lambda = seq(0.01, 4, length.out = 100),
    logpsi = seq(-6, 6, length.out = 100),
    i = c(1, 5)
  ) |>
  mutate(
    ll = future_pmap_dbl(
      .l = list(lambda = lambda, logpsi = logpsi, i = i),
      .f = loglik,
      .progress = TRUE),
    i = factor(i, labels = c("High reliability", "Low reliability"))
  ) 

# Plot of the profile log-likelihood surface
ggplot(plot_df, aes(lambda, logpsi)) +
  geom_raster(aes(fill = ll)) +
  geom_contour(aes(z = ll), col = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red3") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red3") +
  scale_fill_viridis_c() +
  # scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = expression("Profile log-likelihood surface over ("*lambda[i]*", "*psi*")"),
    x = expression(lambda[i]),
    y = expression(log(psi)),
    fill = "Log-likelihood"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_grid(. ~ i)

# Slices of the profile likelihood at fixed psi and lambda
lambda0 <- 
  plot_df |> 
  distinct(lambda) |> 
  slice_min(abs(lambda - 1), with_ties = FALSE) |> 
  pull()
psi0 <- 
  plot_df |> 
  distinct(logpsi) |> 
  slice_min(abs(logpsi - 0), with_ties = FALSE) |> 
  pull()

bind_rows(
  plot_df |>
    filter(logpsi == psi0, lambda < 3) |>
    transmute(i, name = "lambda", value = lambda, ll),
  plot_df |>
    filter(lambda == lambda0, abs(logpsi) < 2.5) |>
    transmute(i, name = "logpsi", value = logpsi, ll)
) |>
  ggplot(aes(value, ll)) +
  geom_line() +
  facet_grid(i ~ name, scales = "free") +
  labs(
    x = NULL, 
    y = "Log-likelihood"
  ) +
  theme_minimal()






