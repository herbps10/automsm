#' Simulate example data for the treatment effect modification NP-MSM.
#'
#' The true conditonal average treatment effect (CATE) is a nonlinear
#' function of \code{X2}, namely \eqn{\psi(X) = \sin(2 \pi X_2)}. This allows
#' the linear working model to be interpreted as the best linear (squared-error)
#' summary of the CATE, while a more flexible working model can
#' track the nonlinear shape.
#'
#' @param N Sample size.
#' @param sigma Conditional standard deviation of the outcome.
#' @param seed Random number seed (optional).
#' @return data frame of simulated data with columns \code{X1}, \code{X2},
#'   \code{A}, and \code{Y}.
#' @importFrom stats runif rbinom rnorm plogis
#' @export
simulate_treatment_effect_modification <- function(
    N = 1e3,
    sigma = 0.1,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  X1 <- stats::runif(N)
  X2 <- stats::runif(N)
  A <- stats::rbinom(N, 1, stats::plogis(X1))

  cate <- sin(2 * pi * X2)

  mu0 <- 0.5 * X1
  mu1 <- 0.5 * X1 + cate
  mu <- ifelse(A == 1, mu1, mu0)
  Y <- stats::rnorm(N, mu, sigma)

  data.frame(X1, X2, A, Y)
}

#' Simulate example data for the conditional dose-response curve NP-MSM.
#' @param N sample size
#' @param treatments number of discrete treatments
#' @param sigma conditional standard deviation of outcome
#' @param seed random number seed (optional)
#' @return data frame of simulated data
#' @importFrom stats runif rnorm
#' @export
simulate_categorical_dose_response <- function(
  N = 1e3,
  treatments = 25,
  sigma = 0.1,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  X1 <- stats::runif(N)
  X2 <- stats::runif(N)
  A <- sample(1:treatments, N, replace = TRUE)
  Y <- stats::rnorm(N, 2 / A, sigma)

  data.frame(X1, X2, A, Y)
}

#' Simulate example data for the longitudinal treatment NP-MSM
#'
#' Simulates a longitudinal data structure with time-varying binary treatments
#' (A_1, ..., A_tau), a single time-varying covariate per time point
#' (L_1, ..., L_tau), and an end-of-study outcome Y that depends
#' linearly on the cumulative treatment duration sum(A_1, ..., A_tau).
#' At each time point, the covariate L_t depends on the previous covariate
#' and treatment, and the treatment A_t is assigned via a logistic model on the
#' history.
#'
#' @param N sample size
#' @param tau number of time pionts
#' @param beta0 intercept of the outcome model
#' @param beta1 coefficient on the cumulative treatment duration in
#'   the outcome model (logit scale)
#' @param seed random number seed (optional)
#' @return data frame of simulated dta with columns L_1, A_1, ..., L_tau, A_tau, Y
#' @importFrom stats runif rbinom rnorm plogis
#' @export
simulate_longitudinal_treatment <- function(
  N = 1e3,
  tau = 3,
  beta0 = -1,
  beta1 = 0.5,
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
  mu <- stats::plogis(beta0 + beta1 * cumtrt + 0.2 * L[, tau])
  Y <- stats::rbinom(
    N,
    1,
    mu
  )

  out <- data.frame(matrix(nrow = N, ncol = 2 * tau + 1))
  colnames(out) <- c(paste0("L", 1:tau), paste0("A", 1:tau), "Y")
  out[, paste0("L", 1:tau)] <- L
  out[, paste0("A", 1:tau)] <- A[, 1:tau]
  out$Y <- Y

  out
}
