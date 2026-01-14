#' Simulate example data for the treatment effect modification NP-MSM.
#' @param N sample size
#' @param sigma conditional standard deviation of outcome
#' @param seed random number seed (optional)
#' @return data frame of simulated data
#' @importFrom stats runif rbinom rnorm plogis
#' @export
simulate_treatment_effect_modification <- function(N = 1e3, sigma = 0.1, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)

  X1 <- stats::runif(N)
  X2 <- stats::runif(N)
  A  <- stats::rbinom(N, 1, stats::plogis(X1))
  mu0 <- 0.5 * X1
  mu1 <- 0.5 * X1 + X2
  mu <- ifelse(A == 1, mu1, mu0)
  Y  <- stats::rnorm(N, mu, sigma)

  data.frame(X1, X2, A, Y)
}
