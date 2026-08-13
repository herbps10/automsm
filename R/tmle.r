#' @importFrom torch torch_tensor optim_lbfgs
#' @noRd
tmle_mle <- function(p, fluctuation_model, ...) {
  epsilon <- torch::torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(epsilon)
  for (iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      target <- fluctuation_model(epsilon, ...)
      #cat(glue::glue("TMLE iteration {iter} target: {as.numeric(target)}\n\n"))
      target$backward(retain_graph = TRUE)
      target
    })
  }
  optimizer$zero_grad()
  epsilon
}

#' Clever covariate for marginal distribution of covariates
#' @noRd
calculate_K <- function(Lm, psi, beta, design_matrix, Minv) {
  p <- dim(Minv)[1]
  Kdot <- batched_dL_dbeta(Lm, psi, beta, design_matrix, p)
  Kdot$matmul(Minv)
}

Q_fluctuation <- function(epsilon, K, Q) {
  Qn <- exp(K$matmul(epsilon) * Q)$sum()
  exp(K$matmul(epsilon) * Q) / Qn
}

dQ_fluctuation_depsilon <- function(epsilon, K, Q) {
  n <- dim(Q)
  Qn <- exp(K$matmul(epsilon) * Q)$sum()
  Q_fluctuation(epsilon, K, Q)$reshape(c(n, 1))$mul(
    K - Q$reshape(c(n, 1))$mul(K)$mul(exp(K * epsilon)) / Qn
  )
}
