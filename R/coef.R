#' Extract coefficients from a "automsm" object
#'
#' @param object An object of class \code{"automsm"}
#' @param estimator Which estimator's influence function to use:
#'   \code{"onestep"} or \code{"tmle"} (default)
#' @param ... Additional arguments (not currently used)
#' @return A named numeric vector of working-model coefficient estimates.
#' @importFrom stats var
#' @export
coef.automsm <- function(object, estimator = c("tmle", "onestep"), ...) {
  estimator <- match.arg(estimator)
  est <- object[[estimator]]
  if(is.null(est) || is.null(est$eif)) {
    stop("No efficient influence function available for estimator '", estimator, "'.", call. = FALSE)
  }

  out <- as.numeric(est$est)
  if(!is.null(object$terms) && length(object$terms) == length(out)) {
    names(out) <- object$terms
  }
  out
}
