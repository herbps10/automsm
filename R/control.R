#' Control parameters for the targeted (TMLE) estimator
#'
#' Constructs the settings governing the iterative targeting step.
#' Pass the result to the \code{tmle} argument of [cate()], [dose_response()],
#' or [longitudinal_dose_response()]. As a shorthand, those arguments also
#' accept a single logical: \code{tmle = FALSE} is equivalent to
#' \code{tmle = tmle_control(enabled = FALSE)}.
#'
#' @param enabled Logical; whether to run the TMLE estimator.
#' @param maxiter Maximum number of outer TMLE iterations.
#' @param fluctuation whether to use a linear TMLE fluctuation model (\code{"linear"}
#'   a logistic one (\code{"logistic"}), or pick based on the outcome type (\code{"auto"}; continuous outcomes default to
#'   linear, and binary outcomes logistic). A logistic fluctuation
#'   model respects the \eqn{[0, 1]} bounds of a binomial outcome regression and
#'   is generally preferable when \code{outcome_type = "binomial"}.
#' @param solver Fluctuation solver: \code{"newton"} (default) uses an
#'   analytic Newton/Iteratively Reweighted Least Squares step, since the fluctuation objective
#'   is a GLM in \eqn{\epsilon} with a closed-form gradient and Hessian. Alternatively,
#'   \code{"lbfgs"} uses L-BFGS with a strong-Wolfe line search over the automatically differentiated
#'   objective. The Newton solver falls back to L-BFGS automatically if it fails or does not decrease
#'   the objective.
#' @param criterion Stopping rule: \code{"eif"} (the default) stops when the EIF estimating equation is approximately solved,
#'   \code{"epsilon"} stops on \eqn{\|\epsilon\|_\infty < \code{tol}}.
#' @param eif_tol Threshold for the \code{"eif"} criterion, in standard-error units.
#' @param tol Threshold for the \code{"epsilon"} criterion. Measured on the \eqn{M^{-1}} scale,
#'   therefore not comparable across problems (this is why \code{"eif"}) is the default).
#' @param eps_max Largest permitted \eqn{\|epsilon\|_\infty} in a single step.
#'   Steps exceeding it are scaled back and a warning is issued, which surfaces a
#'   separated or near-degenerate fluctuation submodel.
#' @param clamp Bound keeping the fluctuated conditional means in
#'   \eqn{[\code{clamp}, 1 - \code{clamp}]} under a logistic fluctuation.
#' @param ridge Relative ridge added to the fluctuation Hessian for numerical stability.
#' @param obj_tol Relative slack when asserting the solver decreased the objective.
#'
#' @return A list of class \code{"automsm_tmle_control"}.
#' @seealso [bayes_control()], [onestep_control()]
#' @family automsm control
#' @export
tmle_control <- function(
  enabled = TRUE,
  maxiter = 25L,
  fluctuation = c("auto", "linear", "logistic"),
  solver = c("newton", "lbfgs"),
  criterion = c("eif", "epsilon"),
  eif_tol = 1e-2,
  tol = NULL,
  eps_max = 25,
  clamp = 1e-3,
  ridge = 1e-10,
  obj_tol = 1e-6
) {
  fluctuation <- match.arg(fluctuation)
  solver <- match.arg(solver)
  criterion <- match.arg(criterion)
  checkmate::assert_flag(enabled)
  checkmate::assert_count(maxiter, positive = TRUE)
  checkmate::assert_number(eif_tol, lower = 0, finite = TRUE)
  checkmate::assert_number(tol, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(eps_max, lower = 0)
  checkmate::assert_number(clamp, lower = 0, upper = 0.5, null.ok = TRUE)
  checkmate::assert_number(ridge, lower = 0)
  checkmate::assert_number(obj_tol, lower = 0)

  structure(list(
    enabled = enabled,
    maxiter = as.integer(maxiter),
    fluctuation = fluctuation,
    solver = solver,
    criterion = criterion,
    eif_tol = eif_tol,
    tol = tol,
    eps_max = eps_max,
    clamp = clamp,
    ridge = ridge,
    obj_tol = obj_tol
  ), class = "automsm_tmle_control")
}


#' Control parameters for the generalied Bayesian TMLE estimator
#'
#' Constructs the settings governing the iterative targeting step,
#' in which a prior is placed on the fluctuation paramter \eqn{epsilon}
#' and combined with the loss function to yield a generalized posterior
#' for the working-model coefficients \eqn{\beta}.
#' Pass the result to the \code{bayes} argument of [cate()], [dose_response()],
#' or [longitudinal_dose_response()]. As a shorthand, those arguments also
#' accept a single logical: \code{bayes = FALSE} is equivalent to
#' \code{bayes = bayes_control(enabled = FALSE)}.
#'
#' The generalized Bayesian estimator is available for the single-time-point estimands only.
#' It requires \code{tmle} to be enabled.
#'
#' @param enabled Logical; whether to run the generalized Bayesian TMLE estimator.
#' @param draws Number of MCMC draws retained per chain.
#' @param warmup Number of initial adaption draws.
#' @param chains Number of MCMC chains.
#' @param prior A function taking the numeric vector of working-model coefficients \eqn{\beta}
#'   (of length \code{p}) and returning a scalar log-density. Defaults to independent standard normals.
#' @param scale Initial proposal covariance for the adaptive sampler. If \code{NULL} (the default),
#'   it is derived from the estimated efficient influence function. Supply a scalar or
#'   length-\eqn{p} vector of variances, or a \eqn{p \times p} covariance matrix to override the default.
#'   Note that this is only for the initial proposal; the sampler adapts toward \code{acc_rate}
#'   during warmup.
#' @param acc_rate Target acceptance rate for the adaptive sampler.
#' @param seed Integer; optional random seed for MCMC algorithm
#' @param labels How to label posterior parameters: \code{"index"} to label by index, or \code{"terms"} to
#'   use the terms from the working model formula
#' @return A list of class \code{"automsm_bayes_control"}.
#' @seealso [tmle_control()], [onestep_control()]
#' @family automsm control
#' @export
bayes_control <- function(
  enabled = TRUE,
  draws = 1e3,
  warmup = 1e3,
  chains = 4L,
  prior = function(beta) sum(stats::dnorm(as.numeric(beta), 0, 1, log = TRUE)),
  scale = NULL,
  acc_rate = 0.3,
  seed = NULL,
  labels = c("index", "terms")
) {
  labels <- match.arg(labels)
  checkmate::assert_flag(enabled)
  checkmate::assert_count(draws, positive = TRUE)
  checkmate::assert_count(warmup, positive = TRUE)
  checkmate::assert_count(chains, positive = TRUE)
  checkmate::assert_function(prior)
  checkmate::assert_numeric(scale, lower = 0, any.missing = FALSE, min.len = 1L, null.ok = TRUE)
  checkmate::assert_number(acc_rate ,lower = 0, upper = 1)
  checkmate::assert_choice(labels, choices = c("index", "terms"))
  structure(list(
    enabled = enabled,
    draws = draws,
    warmup = warmup,
    chains = as.integer(chains),
    prior = prior,
    scale = scale,
    acc_rate = acc_rate
  ), class = "automsm_bayes_control")
}

#' Control parameters for the one-step estimator
#'
#' The one-step estimator corrects the plug-in estimate by the empirical
#' mean of the estimated efficient influence function. It is always computed.
#'
#' @param joint_draws Number of draws from the multivariate normal
#'   approximation to the joint sampling distribution of the one-step estimator.
#' @param seed Option integer seed for the joint draws.
#'
#' @return A list of class \code{"automsm_onestep_control"}.
#' @seealso [tmle_control()], [bayes_control()]
#' @family automsm control
#' @export
onestep_control <- function(
  joint_draws = 1e3,
  seed = NULL
) {
  checkmate::assert_count(joint_draws)
  checkmate::assert_int(seed, null.ok = TRUE)
  structure(list(joint_draws = joint_draws, seed = seed), class = "automsm_onestep_control")
}

#' @noRd
as_control <- function(x, control, cls, arg) {
  if(inherits(x, cls)) return(x)
  if(is.logical(x) && length(x) == 1L && !is.na(x)) return(control(enabled = x))
  stop("`", arg, "` must be TRUE, FALSE, or a `", sub("automsm_", "", sub("_control$", "", cls)), "_control()` object.", call. = FALSE)
}

#' @noRd
as_tmle_control <- function(x) {
  as_control(x, tmle_control, "automsm_tmle_control", "tmle")
}

#' @noRd
as_bayes_control <- function(x) {
  as_control(x, bayes_control, "automsm_bayes_control", "bayes")
}

#' @noRd
as_onestep_control <- function(x) {
  as_control(x, onestep_control, "automsm_onestep_control", "onestep")
}


# ---- Nuisance -----

#' Control parameters for nuisance parameter estimation
#'
#' Constructs the settings governing estimation of the nuisance parameters.
#' Pass the result to the \code{nuisance} argument of [cate()], [dose_response()], or [longitudinal_dose_response()].
#'
#' To supply nuisance estimates computed elsewhere and bypass estimation entirely,
#' use the separate \code{nuisance_estimates} argument instead.
#'
#' @param learners_trt A character vector of \pkg{SuperLearner} libraries for estimating
#'   the propensity scores. For longitudinal estimands the same library is used at every time point.
#' @param learners_outcome A character vector of \pkg{SuperLearner} libraries for
#'   estimating the outcome regressions. For longitudinal estimands the same
#'   library is used for every sequential regression.
#' @param outer_folds Number of folds in the outer cross-fitting loop.
#' @param inner_folds Number of folds for the inner \pkg{SuperLearner} cross-validation
#'   within each outer cross-fitting fold.
#' @param epsilon Truncation bound. Estimated propensity scores are constrained to
#'   \eqn{[\epsilon, 1 - \epsilon]}, and for \code{outcome_type = "binomial"} so are the
#'   estimated conditional means.
#' @param folds Optional pre-constructed fold assignment, as returned by [origami::make_folds()].
#' @param seed Optional integer seed applied to the whole nuisance estimation step.
#'
#' @return A list of class \code{"automsm_nuisance_control"}.
#'
#' @examples
#' # Defaults: 5-fold cross-fitting, main-terms GLM for both nuisances
#' nuisance_control()
#'
#' # Ensemble with reproducible fold assignment
#' nuisance_control(
#'   learners_trt = c("SL.glm", "Sl.ranger", "SL.mean"),
#'   learners_outcome = c("SL.glm", "SL.ranger", "SL.mean"),
#'   outer_folds = 10,
#'   seed = 10016
#' )
#'
#' @seealso [tmle_control()], [bayes_control()], [onestep_control()];
#'   [SuperLearner::SuperLearner()] for the learner libraries;
#'   [origami::make_folds()] for the \code{folds} format.
#'
#' @family automsm control
#' @export
nuisance_control <- function(
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm",
  outer_folds = 5L,
  inner_folds = 5L,
  epsilon = 1e-5,
  folds = NULL,
  seed = NULL
) {
  assert_sl_library(learners_trt, "learners_trt")
  assert_sl_library(learners_outcome, "learners_trt")
  checkmate::assert_count(outer_folds, positive = TRUE)
  checkmate::assert_count(inner_folds, positive = TRUE)
  checkmate::assert_number(epsilon, lower = 0, upper = 0.5)
  checkmate::assert_list(folds, null.ok = TRUE, min.len = 1L)
  checkmate::assert_int(seed, null.ok = TRUE)

  if(inner_folds < 2L) stop("`inner_folds` must be at least 2; SuperLearner requires at least two folds to cross-validate ensemble weights.", call. = FALSE)

  if(outer_folds == 1L && is.null(folds)) {
    warning("`outer_folds = 1` disables cross-fitting: nuisance models are evaluated on the same observations used to fit them.")
  }

  if(is.null(folds) && !is.null(seed)) {
    warning("`seed` has no effect on fold construction when `folds` is supplied; it still seeds SuperLearner's internal cross-validation.")
  }

  if(!is.null(folds)) {
    ok <- all(vapply(folds, function(f) {
      is.list(f) && all(c("training_set", "validation_set") %in% names(f))
    }))
    if(!ok) stop("`folds` must be a list of folds returned by origami::make_folds()", call. = FALSE)
  }

  structure(
    list(learners_trt = learners_trt,
         learners_outcome = learners_outcome,
         outer_folds = as.integer(outer_folds),
         inner_folds = as.integer(inner_folds),
         epsilon = epsilon,
         folds = folds,
         seed = seed),
    class = "automsm_nuisance_control"
  )
}

as_nuisance_control <- function(x) {
  if(inherits(x, "automsm_nuisance_control")) return(x)
  stop("`nuisance` must be a `nuisance_control()` object. To supply pre-computed nuisance estimates, use `nuisance_estimates` instead.", call. = FALSE)
}
