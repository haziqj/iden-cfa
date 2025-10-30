library(tidyverse)
library(lavaan)

## ----- Helper functions to calculate the condition number --------------------
get_cond_number <- function(fit, param = c("anchor", "effcode", "sphere")) {
  param <- match.arg(param)
  eps <- sqrt(.Machine$double.eps)
  
  # Get stuff from lavaan fit
  H  <- lavInspect(fit, "hessian")
  H  <- 0.5 * (H + t(H))
  
  pt <- parTable(fit)
  free_idx <- which(pt$free > 0)
  labs <- pt$label[free_idx]
  
  if (param != "anchor") {
    C <- t(lavInspect(fit, "est")$lambda)
    if (param == "effcode") {
      C[C != 0] <- 1  # the contstraint Jacobian of \sum lambda - p = 0
    }
    if (param == "sphere") {
      C <- 2 * C  # the contstraint Jacobian of \sum lambda^2 - 1 = 0
    } 
    N <- pracma::nullspace(C)
    ix <- which(pt$free > 0 & pt$op == "=~" & pt$rhs %in% colnames(C))
    
    # Full-space projection removing the radial directions from the Lambda block
    P <- diag(ncol(H))
    if (ncol(N) > 0) P[ix, ix] <- N %*% t(N)
    else warning("Degenerate case: constraint Jacobian full rank; tangent space empty.")
    H <- P %*% H %*% P
    H <- 0.5 * (H + t(H))
  }
  
  # Scale-invariant conditioning
  d <- sqrt(pmax(eps, abs(diag(H))))
  S <- H / (d %o% d)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > eps * max(ev)]  # drop numerical zeros (e.g., from projection)
  max(ev) / min(ev)
}

## ----- Compare models --------------------------------------------------------

models <- list(
  anchor1 = "
    visual  =~ 1*x1 + x2 + x3
    textual =~ 1*x4 + x5 + x6
    speed   =~ 1*x7 + x8 + x9
  ",
  anchor2 = "
    visual  =~ NA*x1 + 1*x2 + x3
    textual =~ NA*x4 + 1*x5 + x6
    speed   =~ NA*x7 + 1*x8 + x9
  ",
  anchor3 = "
    visual  =~ NA*x1 + x2 + 1*x3
    textual =~ NA*x4 + x5 + 1*x6
    speed   =~ NA*x7 + x8 + 1*x9
  ",
  fixpsi = "
    visual  =~ NA*x1 + x2 + x3
    textual =~ NA*x4 + x5 + x6
    speed   =~ NA*x7 + x8 + x9
    visual  ~~ 1*visual
    textual ~~ 1*textual
    speed   ~~ 1*speed
  ",
  effcode = "
    visual  =~ NA*x1 + a*x1 + b*x2 + c*x3
    textual =~ NA*x4 + d*x4 + e*x5 + f*x6
    speed   =~ NA*x7 + g*x7 + h*x8 + i*x9
    a + b + c == 3
    d + e + f == 3
    g + h + i == 3",
  sphere = "
    visual  =~ NA*x1 + a*x1 + b*x2 + c*x3
    textual =~ NA*x4 + d*x4 + e*x5 + f*x6
    speed   =~ NA*x7 + g*x7 + h*x8 + i*x9
    n1 := a*a + b*b + c*c
    n2 := d*d + e*e + f*f
    n3 := g*g + h*h + i*i
    n1 == 1
    n2 == 1
    n3 == 1
    a > 0
    d > 0
    g > 0
  "
)

tibble(
  model = models,
  param = c(rep("anchor", 4), "effcode", "sphere")
) |>
  mutate(kappa = map2_dbl(model, param, \(x, y) {
    fit <- cfa(x, HolzingerSwineford1939)
    get_cond_number(fit, y)
  }))
# # A tibble: 6 × 3
#   model        param   kappa
#   <named list> <chr>   <dbl>
# 1 <chr [1]>    anchor   25.5
# 2 <chr [1]>    anchor   60.2
# 3 <chr [1]>    anchor   46.0
# 4 <chr [1]>    anchor   13.4
# 5 <chr [1]>    effcode  17.0
# 6 <chr [1]>    sphere   16.9