#' Tidy method for objects of class "automsm"
#' @param x An object of class \code{"automsm"}
#' @param ... Additional arguments (not currently used)
#' @return A tibble with one row per (estimator, term) combination and columns
#'   \code{estimator}, \code{term}, \code{estimate}, \code{conf.low},
#'   \code{conf.high}, and \code{std.error}
#' @importFrom generics tidy
#' @importFrom stats median quantile sd
#' @importFrom tibble tibble
#' @importFrom purrr map compact list_rbind
#' @exportS3Method tidy automsm
tidy.automsm <- function(x, ...) {
  terms <- x$terms

  tidy_estimator <- function(est, label) {
    if(is.null(est)) return(NULL)

    tibble::tibble(
      estimator = label,
      term = terms,
      estimate = as.numeric(est$est),
      conf.low = as.numeric(est$lower),
      conf.high = as.numeric(est$upper),
      std.error = as.numeric(est$se)
    )
  }

  tidy_samples <- function(draws, label, margin = 3) {
    if(is.null(draws)) return(NULL)
    tibble::tibble(
      estimator = label,
      term = terms,
      estimate = as.numeric(apply(draws, margin, stats::median)),
      conf.low = as.numeric(apply(draws, margin, stats::quantile, 0.025)),
      conf.high = as.numeric(apply(draws, margin, stats::quantile, 0.975)),
      std.error = as.numeric(apply(draws, margin, stats::sd)),
    )
  }

  purrr::list_rbind(purrr::compact(list(
    tidy_estimator(x$onestep, "onestep"),
    tidy_estimator(x$tmle, "tmle"),
    tidy_samples(x$bayes$draws, "bayes")
  )))
}

#' @importFrom generics tidy
#' @export
generics::tidy
