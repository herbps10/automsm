#' Linear working model.
#' @param beta vector of NP-MSM coefficients (must be of type \code{torch::tensor})
#' @param X design matrix (must be of type \code{torch::tensor})
#' @export
working_model_linear <- function(beta, X) {
  X$matmul(beta)
}
