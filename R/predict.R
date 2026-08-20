#' @noRd
AUTOMSM_PREDICT_ESTIMANDS <- c("cate", "dose_response", "longitudinal_dose_response")

#' Predict from a "automsm" object
#'
#' Evaluates the fitted working model \eqn{m_\beta} at supplied design points.
#'
#' @param object An object of class \code{"automsm"}.
#' @param newdata A data frame containing the variables referenced by the working-model formula.
#'   For \code{longitudinal_dose_response} fits that used \code{summary_measures}, this
#'   includes the regime-level summary columns.
#' @param estimator Which estimator to use for predictions:
#'   \code{"tmle"} (default), \code{"onestep"}, or \code{"bayes"}.
#' @param type \code{"point"} for a vector of point predictions;
#'   \code{"joint"} for an \code{ndraws} by \code{nrow(newdata)} matrix of draws;
#'   \code{"rvar"} for a \pkg{posterior} \code{rvar} of length \code{nrow(ndata)}.
#' @param ndraws Number of draws for \code{type="joint"} or \code{"rvar"}. For
#'   \code{estimator = "bayes"}, this subsamples the posterior and is ignored if it exceeds
#'   the number of available draws.
#' @param seed Optional integer seed for the frequentist normal-approximation draws.
#'   Ignored for \code{estimtaor = "bayes"}, whose draws are already fixed by the fit.
#' @param ... Additional arguments (not currently used).
#'
#' @return A numeric vector (\code{"point"}), a numeric matrix with draws in
#'   rows (\code{"joint"}), or a \code{posterior::rvar} (\code{"rvar"}).
#'
#' @seealso [automsm], [posterior::rvar()]
#'
#' @exportS3Method
predict.automsm <- function(
  object,
  newdata,
  estimator = c("tmle", "onestep", "bayes"),
  type = c("point", "joint", "rvar"),
  ndraws = NULL, seed = NULL, ...
) {
  estimator <- match.arg(estimator)
  type <- match.arg(type)
  checkmate::assert_count(ndraws, positive = TRUE, null.ok = TRUE)
  checkmate::assert_int(seed, null.ok = TRUE)

  if(!(object$estimand %in% AUTOMSM_PREDICT_ESTIMANDS)) {
    stop(
      "predict() is not supported for estimand '", object$estimand, "'.",
      call. = FALSE
    )
  }

  if(identical(object$estimand, "longitudinal_dose_response") && identical(estimator, "bayes")) {
    stop("The 'bayes' estimator is not available for the ",
         "'longitudinal_dose_response' estimand.",
         call. = FALSE)
  }

  design <- predict_design(object, newdata)

  beta <- if(identical(estimator, "bayes")) {
    if(is.null(ndraws)) ndraws <- dim(object$bayes_tmle$draws)[1] * dim(object$bayes_tmle$draws)[2]
    predict_beta_bayes(object, if(identical(type, "point")) NULL else ndraws)
  }
  else if(identical(type, "point")) {
    as.numeric(object[[estimator]]$est %||% stop("Estimator '", estimator, "' is not available in this object.", call. = FALSE))
  }
  else {
    if(is.null(ndraws)) ndraws <- 1e3
    predict_beta_freq(object, estimator, ndraws, seed)
  }

  # Bayes + point: summarize the predictions, not beta
  if(identical(estimator, "bayes") && identical(type, "point")) {
    return(apply(predict_eval(object, design, beta), 2, stats::median))
  }

  if(identical(type, "point")) {
    if(anyNA(beta)) {
      stop("Estimator '", estimator, "' contains missing values; cannot predict.", call. = FALSE)
    }
    return(predict_eval(object, design, beta))
  }

  out <- predict_eval(object, design, beta)
  if(identical(type, "rvar")) posterior::rvar(out) else out

}

#' Build the (n_new, 1, d) working-model design from newdata
#'
#' The working model is called with a rank-3 design tensor during fitting, so
#' it must be called with rank 3 here too.
#' @importFrom stats delete.response terms na.pass
#' @noRd
predict_design <- function(fit, newdata) {
  checkmate::assert_data_frame(newdata, min.rows = 1L)

  tf <- stats::delete.response(stats::terms(fit$formula))

  mat <- tryCatch(
    stats::model.matrix(tf, data = newdata, na.action = stats::na.pass),
    error = function(e) {
      stop("Could not build the working-model design from `newdata`: ",
           conditionMessage(e),
           if(identical(fit$estimand, "longitudinal_dose_response")) {
             paste0("\n Note: for the longitudinal estimand the working-model formula may reference the regime-level summary measures ",
                    "(from `summary_measures`), which must be supplied as columns of `newdata`.")
           } else "", call. = FALSE)
    }
  )

  if(nrow(mat) != nrow(newdata)) {
    stop("`newdata` has missing values in the working-model variables.", call. = FALSE)
  }

  if(!identical(colnames(mat), fit$terms)) {
    stop("The design built from `newdata` does not match the fitted working model. \n",
         " fitted: ", paste(fit$terms, collapse = ", "), "\n",
         " newdata: ", paste(colnames(mat), collapse = ", "), "\n",
         "This usually means a factor in `newdata` has different levels than in ",
         "the origianl data, or a variable is missing.", call. = FALSE)
  }

  if(!all(is.finite(mat))) {
    stop("The design built from `newdata` contains non-finite values.",
         call. = FALSE)
  }

  list(
    mat = mat,
    tensor = torch::torch_tensor(mat)$reshape(c(nrow(mat), 1L, ncol(mat)))
  )
}

#' Is the working model affine in beta?
#' @noRd
predict_is_affine <- function(fit, design) {
  if(!identical(fit$p, as.integer(fit$d))) return(FALSE)
  ok <- function(b) {
    got <- tryCatch(
      as.numeric(fit$working_model(torch::torch_tensor(b), design$tensor)),
      error = function(e) NULL
    )
    !is.null(got) && isTRUE(all.equal(got, as.numeric(design$mat %*% b), tolerance = 1e-5))
  }
  ok(rep(0, fit$p)) && ok(seq_len(fit$p) / fit$p) && ok(-seq_len(fit$p) / fit$p)
}

#' Evaluate the working model at one or many beta vectors
#' @noRd
predict_eval <- function(fit, design, beta) {
  if(is.matrix(beta) && nrow(beta) > 1) {
    stopifnot(ncol(beta) == fit$p)
    if(predict_is_affine(fit, design)) return(t(design$mat %*% t(beta)))
    return(t(vapply(seq_len(nrow(beta)), function(i) predict_eval(fit, design, beta[i, ]), numeric(nrow(design$mat)))))
  }
  stopifnot(length(beta) == fit$p)
  out <- fit$working_model(torch::torch_tensor(as.numeric(beta)), design$tensor)
  as.numeric(out$reshape(nrow(design$mat)))
}

#' Beta draws from the frequentist normal approximation
#' @noRd
predict_beta_freq <- function(fit, estimator, ndraws, seed) {
  blk <- fit[[estimator]]

  if(is.null(blk) || is.null(blk$est) || is.null(blk$eif)) {
    stop("Estimator '", estimator, "' is not available in this object.",
         call. = FALSE)
  }
  if(identical(estimator, "tmle") && isFALSE(blk$converged)) {
    stop("The TMLE did not converge (status: ", blk$status %||% "unknown", "), so its estimates are NA. Use `estimator = \"onestep\"`, or refit.", call. = FALSE)
  }

  est <- as.numeric(blk$est)
  eif <- as.matrix(blk$eif)
  if(anyNA(est) || anyNA(eif)) {
    stop("Estimator '", estimator, "' contains missing values; cannot predict.", call. = FALSE)
  }
  sigma <- stats::var(eif) / nrow(eif)
  ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  if(min(ev) < 1e-8 * max(abs(ev), 1)) {
    stop("The estimated sampling covariance is not positive definite; cannot draw from the joint distribution.", call. = FALSE)
  }

  draw <- function() {
    mvtnorm::rmvnorm(ndraws, mean = est, sigma = sigma, checkSymmetry = FALSE, method = "eigen")
  }

  if(is.null(seed)) draw() else withr::with_seed(seed, draw())
}

#' Beta draws from the generalized posterior
#' @noRd
predict_beta_bayes <- function(fit, ndraws) {
  bt <- fit$bayes_tmle
  if(is.null(bt) || is.null(bt$draws)) {
    stop("This fit contains no generalized Bayesian posterior. Refit with `bayes = TRUE` (or `bayes = bayes_control(...)`).", call. = FALSE)
  }

  m <- as.matrix(posterior::as_draws_matrix(bt$draws))
  if(ncol(m) != fit$p) {
    stop("Posterior has ", ncol(m), " variables but the working model has p = ", fit$p, ".", call. = FALSE)
  }
  if(!is.null(ndraws) && ndraws < nrow(m)) {
    m <- m[round(seq(1, nrow(m), length.out = ndraws)), , drop = FALSE]
  }
  m
}

