#' @export
tidy.targeted_msm <- function(x, ...) {
  results <- tibble::tibble(
    estimator = "onestep",
    term      = x$terms,
    estimate  = x$onestep$est,
    conf.low  = x$onestep$lower,
    conf.high = x$onestep$upper,
    std.error = x$onestep$se
  )

  if(!is.null(x$tmle)) {
    results <- rbind(results, tibble::tibble(
      estimator = "tmle",
      term      = x$terms,
      estimate  = x$tmle$est,
      conf.low  = x$tmle$lower,
      conf.high = x$tmle$upper,
      std.error = x$tmle$se
    ))

    if(!is.null(x$tmle$samples)) {
      results <- rbind(results, tibble::tibble(
        estimator = "bayestmle",
        term      = x$terms,
        estimate  = as.numeric(apply(x$tmle$samples, 2, median)),
        conf.low  = as.numeric(apply(x$tmle$samples, 2, quantile, 0.025)),
        conf.high = as.numeric(apply(x$tmle$samples, 2, quantile, 0.975)),
        std.error = as.numeric(apply(x$tmle$samples, 2, sd))
      ))
    }
  }
  results
}
