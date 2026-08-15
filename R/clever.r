#' Batched per-observation clever-covariate directions
#'
#' Returns an (n, K, p) tensor with
#'   d\[i, k, \] = (Minv %*% grad_dL(Lm, psi\[i, \], beta, X\[i, , \]) )\[, k\]
#'
#' Requires \code{psi} and \code{beta} to be leaves with requires_grad = TRUE
#' @noRd
clever_directions <- function(problem, psi, beta, Minv) {
  n <- problem$n
  K <- problem$K
  p <- problem$p
  blocks <- batched_NablaLdot(problem$Lm_fn, psi, beta, problem$design_matrix, p, K)
  Nabla <- torch::torch_stack(blocks, dim = 1)
  d <- torch::torch_matmul(Minv, Nabla$permute(c(2, 1, 3)))
  d$transpose(2, 3)
}
