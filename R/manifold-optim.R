# --- Packages ------------------------------------------------------------
library(ManifoldOptim)
library(Matrix)
library(numDeriv)
library(lavaan)

set.seed(1)

# --- Simulate dataset (your spec) ----------------------------------------
n <- 1000000
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
  
  f1 ~~ 0*f2
"
dat <- lavaan::simulateData(mod, sample.nobs = n)
Y   <- scale(dat, center = TRUE, scale = FALSE) 
p   <- ncol(Y)          # 10
q   <- 2                # factors
S   <- cov(Y)

# --- Parameterization: Lambda = U diag(r), Phi = I, Psi = diag(exp(eta)) --
# Optimize over x = [ vec(U) length p*q ] ⊕ [ z = log r (length q) ] ⊕ [ eta (length p) ]
len_U   <- p * q
len_z   <- q
len_eta <- p
idx_U   <- 1:len_U
idx_z   <- (len_U + 1):(len_U + len_z)
idx_eta <- (len_U + len_z + 1):(len_U + len_z + len_eta)

tx <- function(x) {
  U   <- matrix(x[idx_U], p, q)      # Stiefel-enforced
  z   <- x[idx_z]                    # free reals
  r   <- exp(z)                      # radii > 0
  eta <- x[idx_eta]                  # free reals
  Psi <- diag(exp(eta), p, p)
  Lambda <- U %*% diag(r, q, q)
  list(U = U, r = r, z = z, eta = eta, Psi = Psi, Lambda = Lambda)
}

# --- Negative log-likelihood based on covariance -------------------------
nll <- function(x) {
  par <- tx(x)
  # Phi = I ⇒ Σ = U diag(r^2) U' + Psi
  Sigma <- with(par, U %*% diag(r ^ 2, q, q) %*% t(U) + Psi)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  Ls <- try(chol(Sigma), silent = TRUE)
  if (inherits(Ls, "try-error")) return(1e12)
  logdet <- 2 * sum(log(diag(Ls)))
  X <- backsolve(Ls, t(backsolve(Ls, S, transpose = TRUE)), transpose = TRUE)
  tr <- sum(diag(X))
  n * (logdet + tr)
  # -1 * sum(mvnfast::dmvn(X = Y, mu = rep(0, p), sigma = Ls, log = TRUE, 
  #          isChol = TRUE, ncores = 1))
}

# --- Manifold: Stiefel(p,q) × Euclidean(q + p) ---------------------------
mani.defn <- get.product.defn(
  # get.stiefel.defn(p, q),                 # U
  get.sphere.defn(p, q),
  get.euclidean.defn(len_z + len_eta, 1)  # [z, eta]
)
mani.params   <- get.manifold.params()
solver.params <- get.solver.params(isconvex = FALSE,
                                   maxiter = 2000,
                                   gradtol = 1e-8,
                                   xtol = 1e-10,
                                   ftol = 1e-10)

mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, nll)

# --- Initialization -------------------------------------------------------
# U0: orthonormal columns (Stiefel)
A <- matrix(rnorm(p * q), p, q)
U0 <- qr.Q(qr(A))[, 1:q, drop = FALSE]

# r0: start from PCA-ish scales of S (rough), but keep positive
ev <- eigen(S, symmetric = TRUE)
r0 <- sqrt(pmax(ev$values[1:q], 1e-3))     # rough energy per factor
z0 <- log(r0)

# eta0: uniqueness logs
eta0 <- log(pmax(diag(S), 1e-6) * 0.5)

x0 <- c(as.numeric(U0), z0, eta0)

# --- Optimize -------------------------------------------------------------
res <- manifold.optim(prob,
                      mani.defn,
                      method = "LRBFGS",
                      mani.params = mani.params,
                      solver.params = solver.params,
                      x0 = x0)

cat("Status:", res$status, "  fval:", res$fval, "\n")

parhat <- tx(res$xopt)
U_hat     <- parhat$U
r_hat     <- parhat$r
Lambda_hat<- parhat$Lambda
Psi_hat   <- parhat$Psi
Sigma_hat <- U_hat %*% diag(r_hat^2, q, q) %*% t(U_hat) + Psi_hat

cat("\nU (Stiefel, columns orthonormal):\n"); print(round(U_hat, 3))
cat("\nr (radii):\n"); print(round(r_hat, 3))
cat("\nLambda = U diag(r):\n"); print(round(Lambda_hat, 3))
cat("\nUniquenesses (diag Psi):\n"); print(round(diag(Psi_hat), 3))

# --- Hessian in Grassmann-style tangent (quotient out U -> U R equivalence) ---
x_star <- res$xopt
H_ambient <- numDeriv::hessian(nll, x_star)

# Orthonormal complement for U
Qfull <- qr.Q(qr(U_hat), complete = TRUE)
U_perp <- Qfull[, (q + 1):p, drop = FALSE]   # p × (p-q)

# Tangent basis for U (Grassmann directions only): Delta = U_perp * E_ij
dU_cols <- (p - q) * q
B <- matrix(0, nrow = length(x_star), ncol = dU_cols + len_z + len_eta)

# Fill U directions
col_id <- 1
for (j in 1:q) {
  for (i in 1:(p - q)) {
    D <- matrix(0, p, q)
    D[, j] <- U_perp[, i]
    B[idx_U, col_id] <- as.numeric(D)
    col_id <- col_id + 1
  }
}

# z (log-radii) directions
B[idx_z, (dU_cols + 1):(dU_cols + len_z)] <- diag(len_z)

# eta directions
B[idx_eta, (dU_cols + len_z + 1):(dU_cols + len_z + len_eta)] <- diag(len_eta)

# Tangent Hessian & covariance
H_tangent <- t(B) %*% H_ambient %*% B
H_tangent <- 0.5 * (H_tangent + t(H_tangent))

# Small regularization (exactly equal radii create residual invariance)
eig <- eigen(H_tangent, symmetric = TRUE)
lam <- pmax(eig$values, 1e-9)
Ue  <- eig$vectors
Sigma_tangent <- Ue %*% diag(1/lam, nrow = length(lam)) %*% t(Ue)

# Stable square-root
L_tangent <- t(Ue %*% diag(1/sqrt(lam), nrow = length(lam)))

cat("\nTangent covariance computed (Stiefel + radii).\n")

# --- Varimax rotation on Lambda, then reparameterize back to U diag(r) ---
# 1) Varimax (orthogonal)
vr <- stats::varimax(Lambda_hat)
Lambda_var <- as.matrix(vr$loadings)   # p x q
R_var <- vr$rotmat                     # q x q orthogonal

cat("\nLoadings (varimax-rotated):\n"); print(round(Lambda_var, 3))

# 2) Polar decomposition: Lambda_var = U_pol * S, where S = (Lambda_var' Lambda_var)^{1/2}
M <- crossprod(Lambda_var)             # q x q
eigM <- eigen(M, symmetric = TRUE)
S_half <- eigM$vectors %*% diag(sqrt(pmax(eigM$values, 0)), q, q) %*% t(eigM$vectors)

# Guard against singular M (rare but possible); add tiny ridge if needed
if (min(eigM$values) < 1e-12) {
  eps <- 1e-9
  S_half <- eigM$vectors %*% diag(sqrt(pmax(eigM$values + eps, 0)), q, q) %*% t(eigM$vectors)
}
U_pol <- Lambda_var %*% solve(S_half)  # p x q, should be orthonormal

# 3) Diagonalize S to get radii and a final orthogonal E
eigS <- eigen(S_half, symmetric = TRUE)
E    <- eigS$vectors                   # q x q orthogonal
r_new <- pmax(eigS$values, 0)          # radii >= 0

# Final parameterization consistent with Lambda = U diag(r)
U_new <- U_pol %*% E                   # still Stiefel
Lambda_param <- U_new %*% diag(r_new, q, q)

cat("\nRe-parameterized U (Stiefel) after varimax:\n"); print(round(U_new, 3))
cat("\nRe-parameterized radii r after varimax:\n"); print(round(r_new, 3))
cat("\nLambda_param = U_new diag(r_new):\n"); print(round(Lambda_param, 3))

# Sanity: Sigma is unchanged (up to numerical noise)
Sigma_var <- Lambda_var %*% t(Lambda_var) + Psi_hat
Sigma_param <- Lambda_param %*% t(Lambda_param) + Psi_hat
cat("\n||Sigma_var - Sigma_param||_F:", 
    round(norm(Sigma_var - Sigma_param, "F"), 6), "\n")



efa(dat, nfactors = q)