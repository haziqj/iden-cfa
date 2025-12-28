library(tidyverse)
library(lavaan)
source("R/optim-spherical.R")
source("R/optim-effcode.R")

# Here we have item reliabilities as follows
# rel(y1) = 0.91  [high]
# rel(y2) = 0.80
# rel(y3) = 0.67  [medium]
# rel(y4) = 0.50
# rel(y5) = 0.33  [low]
mod <- "
  f1 =~ y1 + y2 + y3 + y4 + y5 
  f2 =~ y6 + y7 + y8 + y9 + y10
  
  y1 ~~ 0.10*y1
  y2 ~~ 0.25*y2
  y3 ~~ 0.50*y3
  y4 ~~ 1.00*y4
  y5 ~~ 2.00*y5
  
  y6 ~~ 0.10*y6
  y7 ~~ 0.25*y7
  y8 ~~ 0.50*y8
  y9 ~~ 1.00*y9
  y10 ~~ 2.00*y10
  
  f1 ~~ 1*f1
"

sim_fun <- function(b) {
  set.seed(123)
  datlist <- list(
    n50_cor10   = lavaan::simulateData(paste0(mod, "f1 ~~ 0.1*f2"), sample.nobs = 50),
    n5000_cor10 = lavaan::simulateData(paste0(mod, "f1 ~~ 0.1*f2"), sample.nobs = 5000),
    n50_cor30   = lavaan::simulateData(paste0(mod, "f1 ~~ 0.3*f2"), sample.nobs = 50),
    n5000_cor30 = lavaan::simulateData(paste0(mod, "f1 ~~ 0.3*f2"), sample.nobs = 5000),
    n50_cor50   = lavaan::simulateData(paste0(mod, "f1 ~~ 0.5*f2"), sample.nobs = 50),
    n5000_cor50 = lavaan::simulateData(paste0(mod, "f1 ~~ 0.5*f2"), sample.nobs = 5000),
    n50_cor70   = lavaan::simulateData(paste0(mod, "f1 ~~ 0.7*f2"), sample.nobs = 50),
    n5000_cor70 = lavaan::simulateData(paste0(mod, "f1 ~~ 0.7*f2"), sample.nobs = 5000),
    n50_cor90   = lavaan::simulateData(paste0(mod, "f1 ~~ 0.9*f2"), sample.nobs = 50),
    n5000_cor90 = lavaan::simulateData(paste0(mod, "f1 ~~ 0.9*f2"), sample.nobs = 5000)
  )
  
  models <- list(
    anchor_strong = "
    f1 =~ 1*y1 + y2 + y3 + y4 + y5
    f2 =~ 1*y6 + y7 + y8 + y9 + y10
  ",
    # anchor2 = "
    #   f1 =~ NA*y1 + 1*y2 + y3 + y4 + y5
    #   f2 =~ NA*y6 + 1*y7 + y8 + y9 + y10
    # ",
    anchor_medium = "
    f1 =~ NA*y1 + y2 + 1*y3 + y4 + y5
    f2 =~ NA*y6 + y7 + 1*y8 + y9 + y10
  ",
    # anchor4 = "
    #   f1 =~ NA*y1 + y2 + y3 + 1*y4 + y5
    #   f2 =~ NA*y6 + y7 + y8 + 1*y9 + y10
    # ",
    anchor_weak = "
    f1 =~ NA*y1 + y2 + y3 + y4 + 1*y5
    f2 =~ NA*y6 + y7 + y8 + y9 + 1*y10
  ",
    fixpsi = "
    f1 =~ NA*y1 + y2 + y3 + y4 + y5
    f2 =~ NA*y6 + y7 + y8 + y9 + y10
    f1 ~~ 1*f1
    f2 ~~ 1*f2
  ",
    effcode = "
    f1 =~ NA*y1 + a*y1 + b*y2 + c*y3 + d*y4 + e*y5
    f2 =~ NA*y6 + f*y6 + g*y7 + h*y8 + i*y9 + j*y10
    a + b + c + d + e == 5
    f + g + h + i + j == 5
  ",
    sphere = "
    f1 =~ NA*y1 + a*y1 + b*y2 + c*y3 + d*y4 + e*y5
    f2 =~ NA*y6 + f*y6 + g*y7 + h*y8 + i*y9 + j*y10
    n1 := a*a + b*b + c*c + d*d + e*e
    n2 := f*f + g*g + h*h + i*i + j*j
    n1 == 1
    n2 == 1
    a > 0
    f > 0
  "
  )
  
  res_df <- 
    tibble(
      model = models,
      param = names(models)
    ) |>
    expand_grid(
      tibble(
        dat_name = names(datlist),
        dat = datlist
      )
    ) |>
    mutate(kappa = pmap_dbl(list(model, param, dat), \(x, y, z) {
      fit <- cfa(x, z)
      S <- cov(z)
      n <- nrow(z)
      
      if (y == "effcode") {
        kappa <- planar_cfa(S, n)
      } else if (y == "sphere") {
        kappa <- spherical_cfa(S, n)
      } else {
        eps <- sqrt(.Machine$double.eps)
        H  <- lavInspect(fit, "hessian")
        H  <- 0.5 * (H + t(H))
        d <- sqrt(pmax(eps, abs(diag(H))))
        S <- H / (d %o% d)
        ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        ev <- ev[ev > eps * max(ev)]  # drop numerical zeros (e.g., from projection)
        kappa <- max(ev) / min(ev)
      }
      kappa
    })) 
  
  res_df |>
    select(param, dat_name, kappa) |>
    separate(dat_name, into = c("n", "cor"), sep = "_") |>
    filter(n == "n50") |>
    pivot_wider(names_from = c(cor), values_from = c(kappa))
  
}

library(furrr)
plan(multisession, workers = parallel::detectCores() - 1)
B <- 100  # no. of simulations
res <- future_map(seq_len(B), sim_fun, .progress = TRUE, 
                  .options = furrr_options(seed = TRUE))

bind_rows(res, .id = "i") |> 
  summarise(
    across(starts_with("cor"), mean),
    .by = c(param)
  )
