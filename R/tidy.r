#' Tidy method for objects of class "targeted_msm"
#' @param x An object of class \code{"targeted_msm"}
#' @param ... Additional arguments (not currently used)
#' @return A tibble with one row per (estimator, term) combination and columns
#'   \code{estimator}, \code{term}, \code{estimate}, \code{conf.low},
#'   \code{conf.high}, and \code{std.error}
#' @importFrom generics tidy
#' @importFrom stats median quantile sd
#' @importFrom tibble tibble
#' @importFrom purrr map compact list_rbind
#' @exportS3Method tidy targeted_msm
tidy.targeted_msm <- function(x, ...) {
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

  tidy_samples <- function(samples, label, margin = 3) {
    if(is.null(samples)) return(NULL)
    tibble::tibble(
      estimator = label,
      term = terms,
      estimate = as.numeric(apply(samples, margin, stats::median)),
      conf.low = as.numeric(apply(samples, margin, stats::quantile, 0.025)),
      conf.high = as.numeric(apply(samples, margin, stats::quantile, 0.975)),
      std.error = as.numeric(apply(samples, margin, stats::sd)),
    )
  }

  purrr::list_rbind(purrr::compact(list(
    tidy_estimator(x$onestep, "onestep"),
    tidy_estimator(x$tmle, "tmle"),
    tidy_samples(x$tmle$samples, "bayestmle")
  )))
}

#' @importFrom generics tidy
#' @export
generics::tidy
