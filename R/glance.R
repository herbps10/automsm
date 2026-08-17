#' Glance at a "automsm" object
#'
#' @param x An object of class \code{"automsm"}
#' @param ... Additional arguments (not currently used).
#' @return A one-row tibble summarizing the fitted model
#' @importFrom generics glance
#' @importFrom tibble tibble
#' @export
glance.automsm <- function(x, ...) {
  out <- tibble::tibble(
    estimand = x$estimand,
    n = x$n,
    p = x$p
  )

  if(!is.null(x$tau)) {
    out$tau <- x$tau
  }

  out
}


#' @importFrom generics glance
#' @export
generics::glance
