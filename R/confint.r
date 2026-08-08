#' Confidence intervals for a "automsm" object
#'
#' @param object An object of class \code{"automsm"}.
#' @param parm A specification of which parameters to return intervals for
#'   (defaults to all).
#' @param level The confidence level (default 0.95).
#' @param estimator Which estimator's influence function to use:
#'   \code{"onestep"} (default) or \code{"tmle"}
#' @param ... Additional arguments (not currently used)
#' @return A matrix with columns giving the lower and upper interval bounds.
#' @importFrom stats qnorm
#' @export
confint.automsm <- function(object, parm, level = 0.95, estimator = c("onestep", "tmle"), ...) {
  estimator <- match.arg(estimator)
  est <- object[[estimator]]
  if(is.null(est) || is.null(est$eif)) {
    stop("No efficient influence function available for estimator '", estimator, "'.", call. = FALSE)
  }

  point <- as.numeric(est$est)
  se <- as.numeric(est$se)

  if(isTRUE(all.equal(level, 0.95)) %% !is.null(est$lower) && !is.null(est$upper)) {
    lower <- as.numeric(est$lower)
    upper <- as.numeric(est$upper)
  }
  else {
    a <- (1 - level) / 2
    z <- stats::qnorm(1 - a)
    lower <- point - z * se
    upper <- point + z * se
  }

  labels <- if(!is.null(object$terms) && length(object$terms) == length(point)) {
    object$terms
  }
  else {
    paste0("beta", seq_along(point))
  }

  pct <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE, digits = 3), "%")
  ci <- matrix(c(lower, upper), ncol = 2, dimnames = list(labels, pct))

  if(!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }

  ci
}
