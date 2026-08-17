#' @noRd
beta_variable_names <- function(problem, labels = c("index", "terms")) {
  labels <- match.arg(labels)
  p <- problem$p
  if(labels == "terms") {
    if(length(problem$terms) != p) {
      warning("Cannot label coefficients by design-matrix term: ",
              "the working model has p = ", p, " coefficients but the design has ",
              length(problem$terms), " columns. Using indexed names.",
              call. = FALSE)
      return(sprintf("beta[%d]", seq_len(p)))
    }
    return(problem$terms)
  }
  sprintf("beta[%d]", seq_len(p))
}

#' posterior::draws_array from an MCMC sample
#'
#' @param samples array (chains x iterations x variables) made by run_bayes_tmle
#' @param variables character vector of variable names
#' @noRd
new_draws_array <- function(samples, variables, terms = NULL) {
  stopifnot(length(dim(samples)) == 3L, dim(samples)[3] == length(variables))
  # Change our array layout (chain x iteration x variable) to (iteration x chain x variable)
  x <- aperm(samples, c(2L, 1L, 3L))
  dimnames(x) <- list(
    iteration = as.character(seq_len(dim(x)[1])),
    chain = as.character(seq_len(dim(x)[2])),
    variable = variables
  )
  out <- posterior::as_draws_array(x)
  if(!is.null(terms)) attr(out, "automsm_terms") <- terms
  out
}

#' Convergence diagnostics
#' @importFrom stats quantile
#' @noRd
bayes_diagnostics <- function(draws, warn = TRUE) {
  diag <- posterior::summarise_draws(
    draws, "mean", "median", "sd", "mad",
    q5 = ~stats::quantile(.x, 0.05), q95 = ~stats::quantile(.x, 0.95),
    "rhat", "ess_bulk", "ess_tail"
  )
  if(warn) {
    bad_rhat <- diag$variable[!is.na(diag$rhat) & diag$rhat > 1.01]
    if(length(bad_rhat)) {
      warning("Generalized posterior may not have converged: R-hat > 1.01 for ",
              paste(bad_rhat, collapse = ", "),
              ". Consider increasing `draws`/`warmup`, or tuning `scale`.", call. = FALSE)
    }

    low_ess <- diag$variable[!is.na(diag$ess_bulk) & diag$ess_bulk < 100]
    if(length(low_ess)) {
      warning("Low effective sample size (bulk ESS < 100) for ",
              paste(low_ess, collapse = ", "),
              ". The MCMC random-walk sampler may be mixing if `scale` is poorly matched to the posterior.",
              call. = FALSE)
    }
  }
  diag
}

#' Coerce a fitted NP-MSM to posterior draws
#'
#' Extracts the generalized Bayesian posterior from an \code{automsm} fit produced
#' by [cate()] or [dose_response()] with the generalized Bayesian estimator enabled.
#' Returns the draws in one of the \pkg{posterior} draws formats for use with the
#' wider Bayesian worfkow.
#'
#' @param x An object of class \code{"automsm"}
#' @param parameter Which parameter to extract: \code{"beta"} (default) for
#'   the working model coefficients, or \code{"epsilon"} for the fluctuation parameter
#'   together with the log-density (\code{lp__}).
#' @param ... Not currently used.
#'
#' @details
#' Default variable names are \code{beta[1]}, ..., \code{beta[p]}. Use
#' \code{bayes_control(labels = "terms")} to name coefficients by the design matrix term
#' instead. The term labels are available via the \code{"automsm_terms"} attribute either way.
#'
#' @return A \pkg{posterior} draws object.
#' @seealso [bayes_control()], [automsm]
#' @export
as_draws_array.automsm <- function(x, parameter = c("beta", "epsilon"), ...) {
  parameter <- match.arg(parameter)
  if(is.null(x$bayes_tmle)) {
    stop("This fit contains no generalized Bayesian posterior. Refit with ",
         "`bayes = TRUE` (or `bayes = bayes_control(...)`.", call. = FALSE)
  }
  switch(parameter, beta = x$bayes_tmle$draws, epsilon = x$bayes_tmle$epsilon_draws)
}

#' @rdname as_draws_array.automsm
#' @export
as_draws.automsm <- function(x, ...) as_draws_array.automsm(x, ...)

#' @rdname as_draws_array.automsm
#' @export
as_draws_matrix.automsm <- function(x, ...) {
  posterior::as_draws_matrix(as_draws_array.automsm(x, ...))
}

#' @rdname as_draws_array.automsm
#' @export
as_draws_df.automsm <- function(x, ...) {
  posterior::as_draws_df(as_draws_array.automsm(x, ...))
}

#' @rdname as_draws_array.automsm
#' @export
as_draws_list.automsm <- function(x, ...) {
  posterior::as_draws_list(as_draws_array.automsm(x, ...))
}

#' @rdname as_draws_array.automsm
#' @export
as_draws_rvars.automsm <- function(x, ...) {
  posterior::as_draws_rvars(as_draws_array.automsm(x, ...))
}
