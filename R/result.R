#' Canonical field order for \code{automsm} objects
#'
#' Every estimand returns exactly these fields, in this order.
#' \code{regimes} is NULL for single-time-point estimands.
#' \code{tmle}/\code{bayes_tmle} are \code{NULL} when
#' the corresponding estimator is disabled.
#' Fields are never absent.
#' @noRd

AUTOMSM_FIELDS <- c(
  "estimand", "p", "d", "n", "tau",
  "formula", "working_model", "loss", "terms",
  "nuisance_control",
  "nuisance_estimates",
  "regimes",
  "plugin", "onestep", "tmle", "bayes_tmle"
)

#' @noRd
new_automsm <- function(problem, base, tmle = NULL, bayes_tmle = NULL, nuisance) {
  res <- list(
    estimand           = problem$estimand,
    p                  = problem$p,
    d                  = problem$d,
    n                  = problem$n,
    tau                = problem$tau,
    formula            = problem$formula,
    working_model      = problem$working_model,
    loss               = problem$loss,
    terms              = problem$terms,
    nuisance_control   = nuisance,
    nuisance_estimates = problem$nuisance_estimates,
    regimes            = problem$aux$regimes,
    plugin             = base$plugin,
    onestep            = base$onestep,
    tmle               = tmle,
    bayes_tmle         = bayes_tmle
  )
  stopifnot(identical(names(res), AUTOMSM_FIELDS))
  structure(res, class = "automsm")
}
