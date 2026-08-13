#' Simulate example data for the Conditional Average Treatment Effect NP-MSM.
#'
#' The true conditonal average treatment effect (CATE) is a nonlinear
#' function of \code{X2}, namely \eqn{\psi(X) = \sin(2 \pi X_2)}. This allows
#' the linear working model to be interpreted as the best linear (squared-error)
#' summary of the CATE, while a more flexible working model can
#' track the nonlinear shape.
#'
#' @param N Sample size.
#' @param nonlinear Whether to generate a nonlinear (\code{TRUE}; default) or linear (\code{FALSE}) outcome.
#' @param binary Whether to simulate binary (\code{TRUE}; default) or continuous (\code{FALSE}) outcomes.
#' @param sigma Conditional standard deviation of the outcome (only used if \code{continuous = TRUE}).
#' @param seed Random number seed (optional).
#' @return data frame of simulated data with columns \code{X1}, \code{X2},
#'   \code{A}, and \code{Y}.
#' @importFrom stats runif rbinom rnorm plogis
#' @export
simulate_cate <- function(
    N = 1e3,
    nonlinear = TRUE,
    binary = TRUE,
    sigma = 0.1,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  X1 <- stats::runif(N)
  X2 <- stats::runif(N)
  A <- stats::rbinom(N, 1, stats::plogis(X1))

  if(isTRUE(nonlinear)) {
    cate <- sin(2 * pi * X2)
  }
  else {
    cate <- 0.25 * X2
  }

  mu0 <- 0.5 * X1
  mu1 <- 0.5 * X1 + cate

  if(isTRUE(binary)) {
    mu0 <- plogis(mu0)
    mu1 <- plogis(mu1)
  }

  mu <- ifelse(A == 1, mu1, mu0)

  Y <- if(isTRUE(binary)) stats::rbinom(N, 1, mu) else stats::rnorm(N, mu, sigma)

  data.frame(X1, X2, A, Y)
}

#' Simulate example data for the conditional dose-response curve NP-MSM.
#'
#' @param N sample size
#' @param treatments number of discrete treatments
#' @param nonlinear Whether to generate a nonlinear (\code{TRUE}; default) or linear (\code{FALSE}) outcome.
#' @param binary Whether to simulate binary (\code{TRUE}; default) or continuous (\code{FALSE}) outcomes.
#' @param sigma Conditional standard deviation of the outcome (only used if \code{continuous = TRUE}).
#' @param seed random number seed (optional)
#' @return data frame of simulated data
#' @importFrom stats runif rnorm rbinom
#' @export
simulate_dose_response <- function(
  N = 1e3,
  treatments = 25,
  nonlinear = TRUE,
  binary = TRUE,
  sigma = 0.1,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  X1 <- stats::runif(N)
  X2 <- stats::runif(N)
  A <- sample(1:treatments, N, replace = TRUE)

  mu <- if(isTRUE(nonlinear)) 2 / A else A

  if(isTRUE(binary)) {
    mu <- plogis(mu)
    Y <- stats::rbinom(N, 1, mu)
  }
  else {
    Y <- stats::rnorm(N, mu, sigma)
  }

  data.frame(X1, X2, A, Y)
}

#' Simulate example data for the longitudinal dose-response NP-MSM
#'
#' Simulates a longitudinal data structure with time-varying binary treatments
#' (A_1, ..., A_tau), a single time-varying covariate per time point
#' (L_1, ..., L_tau), and an end-of-study outcome Y that depends
#' linearly on the cumulative treatment duration sum(A_1, ..., A_tau).
#' At each time point, the covariate L_t depends on the previous covariate
#' and treatment, and the treatment A_t is assigned via a logistic model on the
#' history.
#'
#' @param N sample size.
#' @param tau number of time points.
#' @param beta0 intercept of the outcome model.
#' @param beta1 coefficient on the cumulative treatment duration in
#'   the outcome model (logit scale).
#' @param binary whether to simulate binary (default) or continuous outcomes.
#' @param seed random number seed (optional).
#' @return data frame of simulated dta with columns L_1, A_1, ..., L_tau, A_tau, Y
#' @importFrom stats runif rbinom rnorm plogis
#' @export
simulate_longitudinal_dose_response <- function(
  N = 1e3,
  tau = 3,
  beta0 = -1,
  beta1 = 0.5,
  binary = TRUE,
  sigma = 0.1,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  L <- matrix(NA, nrow = N, ncol = tau)
  A <- matrix(NA, nrow = N, ncol = tau)

  L[, 1] <- stats::rnorm(N)
  pscore <- stats::plogis(1 + 0.3 * L[, 1])
  A[, 1] <- stats::rbinom(N, 1, pscore)

  if (tau >= 2) {
    for (t in 2:tau) {
      L[, t] <- stats::rnorm(N, mean = 0.4 * L[, t - 1] + 0.5 * A[, t - 1])
      A[, t] <- stats::rbinom(
        N,
        1,
        stats::plogis(0.3 * L[, t] + 0.4 * A[, t - 1])
      )
    }
  }

  cumtrt <- rowSums(A)
  if(binary == TRUE) {
    mu <- stats::plogis(beta0 + beta1 * cumtrt + 0.2 * L[, tau])
    Y <- stats::rbinom(N, 1, mu)
  }
  else {
    mu <- beta0 + beta1 * cumtrt + 0.2 * L[, tau]
    Y <- stats::rnorm(N, mu, sigma)
  }

  out <- data.frame(matrix(nrow = N, ncol = 2 * tau + 1))
  colnames(out) <- c(paste0("L", 1:tau), paste0("A", 1:tau), "Y")
  out[, paste0("L", 1:tau)] <- L
  out[, paste0("A", 1:tau)] <- A[, 1:tau]
  out$Y <- Y

  out
}
