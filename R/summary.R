#' Summary method for objects of class "automsm"
#' @param object An object of class \code{"automsm"}
#' @param ... Additional arguments (not currently used)
#' @return An object of class \code{"summary.automsm"}: a
#'   list containing the esetimand label, sample size, working-model
#'   term names, and a list of coefficient tables (one per available
#'   estimator: plug-in, one-step, and optionally TMLE).
#' @export
summary.automsm <- function(object, ...) {
  # Build a coefficient table for each available estimator
  estimators <- list()
  if(!is.null(object$onestep)) {
    estimators[["One-step"]] <- msm_coef_table(object$onestep, object$terms)
  }
  if(!is.null(object$tmle)) {
    estimators[["TMLE"]] <- msm_coef_table(object$tmle, object$terms)
  }

  structure(
    list(
      estimand = object$estimand,
      label = automsm_label(object$estimand),
      n = object$n,
      terms = object$terms,
      estimators = estimators
    ),
    class = "summary.automsm"
  )
}

#' Assemble a coefficient table from an estimator sub-list
#'
#' @param est A list with numeric components \code{est}, \code{se},
#'   \code{lower}, \code{upper}.
#' @param terms Optional character vector of coefficient names.
#' @noRd
msm_coef_table <- function(est, terms = NULL) {
  p <- length(est$est)
  if(is.null(terms) || length(terms) != p) {
    terms <- paste0("beta", seq_len(p))
  }

  data.frame(
    term = terms,
    est = as.numeric(est$est),
    se = as.numeric(est$se),
    lower = as.numeric(est$lower),
    upper = as.numeric(est$upper),
    stringsAsFactors = FALSE
  )
}

#' Print method for objects of class "summary.automsm"
#'
#' @param x An object of class \code{"summary.automsm"}.
#' @param digits Number of digits for the coefficient tables.
#' @param ... Additional arguments (not currently used).
#' @return \code{x}, invisibly.
#' @importFrom cli cli_h1 cli_h2 cli_text cli_verbatim
#' @export
print.summary.automsm <- function(x, digits = 3, ...) {
  cli::cli_h1("automsm: {x$label}")
  cli::cli_text("{.field n} = {x$n}")

  if(length(x$estimators) == 0) {
    cli::cli_text("{.emph No fitted estimators available.}")
  }

  for(nm in names(x$estimators)) {
    cli::cli_h2("{nm} estimator")
    tbl <- x$estimators[[nm]]
    print_coef_table(tbl, digits = digits)
  }

  invisible(x)
}

#' Pretty-print a coefficient table
#' @noRd
print_coef_table <- function(tbl, digits = 3) {
  out <- data.frame(
    Term = tbl$term,
    Est = round(tbl$est, digits),
    SE = round(tbl$se, digits),
    `2.5%` = round(tbl$lower, digits),
    `97.5%` = round(tbl$upper, digits),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cli::cli_verbatim(paste(
    utils::capture.output(print(out, row.names = FALSE)),
    collapse = "\n"
  ))
}
