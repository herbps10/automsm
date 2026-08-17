#' Minimise a fluctuation objective over epsilon
#'
#' \code{obj_fn} is a closure of \code{epsilon} only.
#'
#' @importFrom torch torch_tensor optim_lbfgs
#' @noRd
tmle_mle <- function(p, obj_fn) {
  epsilon <- torch::torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(epsilon, max_iter = 20)

  for (iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      target <- obj_fn(epsilon)
      target$backward(retain_graph = TRUE)
      target
    })
  }
  optimizer$zero_grad()
  epsilon
}

#' Clever covariate for marginal distribution of covariates
#' @noRd
calculate_K <- function(Lm, psi, beta, design_matrix, Minv, batched_beta) {
  p <- dim(Minv)[1]
  Kdot <- batched_dL_dbeta(Lm, psi, beta, design_matrix, p, batched_beta)
  Kdot$matmul(Minv)
}

Q_fluctuation <- function(epsilon, K, Q) {
  u <- K$matmul(epsilon)
  u <- u - u$max() # stabilizes estimation
  w <- torch::torch_exp(u)$mul(Q)
  w$divide(w$sum())
}

dQ_fluctuation_depsilon <- function(epsilon, K, Q) {
  n <- dim(Q)
  Qt <- Q_fluctuation(epsilon, K, Q)
  Qt$reshape(c(n, 1))$mul(K - Qt$reshape(c(1, n))$matmul(K))
}
