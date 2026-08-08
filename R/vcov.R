#' Variance-covariance matrix for a "targeted_msm" object
#'
#' Computes the estimated covariance matrix as the empirical variance of the
#' estimated efficient influence function divided by sample size.
#'
#' @param object An object of class \code{"targeted_msm"}.
#' @param estimator Which estimator's influence function to use:
#'   \code{"onestep"} or \code{"tmle"} (default).
#' @param ... Additional arguments (not currently used)
#' @return A p x p covariance matrix.
#' @importFrom stats var
#' @export
vcov.targeted_msm <- function(object, estimator = c("tmle", "onestep"), ...) {
  estimator <- match.arg(estimator)
  est <- object[[estimator]]
  if(is.null(est) || is.null(est$eif)) {
    stop("No efficient influence function available for estimator '", estimator, "'.", call. = FALSE)
  }

  V <- stats::var(as.matrix(est$eif)) / object$n

  labels <- if(!is.null(object$terms) && nrow(V) == length(object$terms)) {
    object$terms
  }
  else {
    paste0("beta", seq_len(nrow(V)))
  }
  dimnames(V) <- list(labels, labels)
  V
}
