# ----- Simulation study 1 DGP -----
sim1_data <- function(n = 200L, seed = 10016) {
  set.seed(seed)
  X <- matrix(rnorm(4 * n), ncol = 4)
  colnames(X) <- paste0("X", 1:4)
  d <- as.data.frame(X)
  d$pi <- plogis(0.5 * d$X1 - 0.5 * d$X2 + 0.2 * d$X3 - 0.1 * d$X4)
  d$A <- rbinom(n, 1L, d$pi)
  d$mu0 <- plogis(0.5 * d$X2 - 0.5 * d$X3)
  d$mu1 <- plogis(0.5 * d$X2 - 0.5 * d$X3 + 1 + 0.5 * d$X4)
  d$mu <- ifelse(d$A == 1, d$mu1, d$mu0)
  d$condvar <- d$mu * (1 - d$mu)
  d$condvar0 <- d$mu0 * (1 - d$mu0)
  d$condvar1 <- d$mu1 * (1 - d$mu1)
  d$Y <- rbinom(n, 1L, d$mu)
  d
}

# DGP with continuous outcomes and different variance in each arm
sim2_data <- function(n = 200L, seed = 10016) {
  set.seed(seed)
  X <- matrix(rnorm(4 * n), ncol = 4)
  colnames(X) <- paste0("X", 1:4)
  d <- as.data.frame(X)
  d$pi <- plogis(0.5 * d$X1 - 0.5 * d$X2 + 0.2 * d$X3 - 0.1 * d$X4)
  d$A <- rbinom(n, 1L, d$pi)
  d$mu0 <- 0.5 * d$X2 - 0.5 * d$X3
  d$mu1 <- 0.5 * d$X2 - 0.5 * d$X3 + 1 + 0.5 * d$X4
  d$mu <- ifelse(d$A == 1, d$mu1, d$mu0)
  d$condvar0 <- 0.1^2
  d$condvar1 <- 0.25^2
  d$condvar <- ifelse(d$A == 1, d$condvar1, d$condvar0)
  d$Y <- rnorm(
    n,
    mean = ifelse(d$A == 1, d$mu1, d$mu0),
    sd = sqrt(ifelse(d$A == 1, d$condvar1, d$condvar0))
  )
  d
}

oracle_nuisance_cate <- function(data) {
  mu <- ifelse(data$A == 1L, data$mu1, data$mu0)
  list(
    pi = data$pi,
    mu = mu,
    mu0 = data$mu0,
    mu1 = data$mu1,
    condvar = data$condvar,
    condvar0 = data$condvar0,
    condvar1 = data$condvar1
  )
}

oracle_nuisance_dose_response <- function(data) {
  mu <- ifelse(data$A == 1L, data$mu1, data$mu0)
  list(
    pi = ifelse(data$A == 1L, data$pi, 1 - data$pi),
    mu = mu,
    mu_a = cbind(data$mu0, data$mu1),
    pi_a = cbind(1 - data$pi, data$pi),
    condvar = data$condvar,
    condvar_a = cbind(data$condvar0, data$condvar1)
  )
}

# ----- Oracle nuisance -----

cate_delta <- function(data, nuisance) {
  H <- ifelse(data$A == 1L, 1 / nuisance$pi, -1 / (1 - nuisance$pi))
  H * (data$Y - nuisance$mu)
}

# Closed-form EIF for the linear working model + squared-error loss
analytic_eif_linear <- function(Xmat, psi, beta, Delta) {
  n <- nrow(Xmat)
  Sigma <- crossprod(Xmat) / n
  r <- Delta + (psi - as.vector(Xmat %*% beta))
  (r * Xmat) %*% solve(Sigma)
}

# ----- Custom expectations -----

expect_solves_eif <- function(eif, frac = 1e-2) {
  eif <- as.matrix(eif)
  n <- nrow(eif)
  se <- apply(eif, 2, stats::sd) / sqrt(n)
  testthat::expect_lt(max(abs(colMeans(eif))), frac * max(se))
}

expect_tensor_close <- function(actual, expected, tolerance = 1e-6) {
  expect_equal(
    as.array(actual$detach()$cpu()),
    as.array(expected),
    tolerance = tolerance
  )
}

#' Jacobian of a torch function via autograd, one row at a time
#' @noRd
autograd_jacobian <- function(f, x, dtype = torch::torch_double()) {
  e <- torch::torch_tensor(x, dtype = dtype, requires_grad = TRUE)
  out <- f(e)
  n <- length(out)
  t(vapply(seq_len(n), function(i) {
    as.numeric(torch::autograd_grad(out[i], e, retain_graph = TRUE)[[1]])
  }, numeric(length(x))))
}

numeric_jacobian <- function(f, x, h = 1e-6) {
  f0 <- f(x)
  J <- array(0, dim = c(length(f0), length(x)))
  for(j in seq_along(x)) {
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

leaf <- function(x) torch::torch_tensor(x)$detach()$clone()$requires_grad_(TRUE)
