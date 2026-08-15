#' Linear working model.
#' @param beta vector of NP-MSM coefficients (must be of type \code{torch::tensor})
#' @param X design matrix (must be of type \code{torch::tensor})
#' @export
working_model_linear <- function(beta, X) {
  if(beta$dim() == 1L) {
    X$matmul(beta)
  } else {
    torch::torch_einsum("nkd,nd->nk", list(X, beta))
  }
}
