

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
