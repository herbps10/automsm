tmle_control <- function(enabled = TRUE, maxiter = 25L, linear = TRUE, tol = NULL) {
  list(enabled = enabled, maxiter = as.integer(maxiter), linear = linear, tol = tol)
}

bayes_control <- function(enabled = FALSE, draws = 1e3, chains = 4L,
                          prior = function(beta) {
                            sum(stats::dnorm(as.numeric(beta), 0, 1, log = TRUE))
                          }, scale = 1e-3, acc_rate = 0.3) {
  list(enabled = enabled, draws = draws, chains = as.integer(chains), prior = prior, scale = scale, acc_rate = acc_rate)
}

onestep_control <- function(joint_draws = 1e3, seed = NULL) {
  list(joint_draws = joint_draws, seed = seed)
}
