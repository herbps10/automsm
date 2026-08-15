#' Shared arguments for automsm estimators
#'
#' @param data A \code{data.frame} containing all referenced columns.
#' @param formula A model formula specifying the marginal structural working-model design matrix.
#' @param outcome_type Outcome type, either \code{"continuous"} or \code{"binomial"}.
#' @param loss The marginal structural model loss function measuring the fidelity of the
#'   working model to the target functional (e.g. \code{loss_squared_error}).
#' @param working_model The marginal structural model working model (e.g. \code{working_model_linear}).
#' @param p Integer giving the dimension of the working-model parameter
#'   \eqn{\beta} (the number of working-model coefficients). If \code{NULL}
#'   (the default), \code{p} is inferred from the number of columns in the
#'   design matrix produced by \code{formula}. When a working model whose
#'   coefficient dimension does not match the \code{formula}-based design, \code{p}
#'   should be supplied explicitly to match the working model's true number of coefficients.
#' @param nuisance Settings for nuisance estimation: a [nuisance_control()] object.
#' @param nuisance_estimates Optional list of pre-computed nuisance parameters. If \code{NULL}
#'   (the default), nuisance parameters are estimated internally via cross-fitting.
#' @param tmle Settings for the targeted estimator: either a logical or a [tmle_control()] object.
#' @param bayes Settings for the generalized Bayesian estimator: either a logical or a [bayes_control()] object.
#' @param onestep Settings for the one-step estimator: a [onestep_control()] object.
#'
#' @name automsm_shared_params
#' @keywords internal
NULL

#' Fitted non-parametric marginal structural model
#'
#' Objects returned by [cate()], [dose_response()],
#' and [longitudinal_dose_response()]. All three return the same fields,
#' in the same order. Fields that do not apply are \code{NULL} rather than absent.
#'
#' @section Components:
#'   \describe{
#'     \item{estimand}{Character string identifying the estimand.}
#'     \item{p}{Number of working-model coefficients.}
#'     \item{d}{Number of columns of the working-model design matrix. Equals
#'       \code{p} unless a working model with a differing coefficnet count was supplied.}
#'     \item{n}{Sample size.}
#'     \item{tau}{Number of treatment timepoints; \code{1} for single-time-point estimands.}
#'     \item{formula, working_model, loss, terms}{The working model as specified, and the design-matrix term names.}
#'     \item{learners_trt, learners_outcome}{The \pkg{SuperLearner} libraries used.}
#'     \item{nuisance}{The estimated or supplied nuisance parameters.}
#'     \item{regimes}{The treatment trajectories over which the model is defined,
#'       or \code{NULL} for the single-time-point estimands.}
#'     \item{plugin}{A list with the plug-in point estimate (\code{est}).}
#'     \item{onestep}{A list with the one-step point estimate (\code{est}), standard
#'       errors (\code{se}), confidence-interval bounds (\code{lower}, \code{upper}),
#'       the estimated efficient influence function (\code{eif}), and draws from the joint sampling distribution
#'       (\code{joint_draws}).}
#'     \item{tmle}{List with the same fields as \code{onestep} plus \code{converged} (logical) and \code{iter} (TMLE iterations used). \code{NULL}
#'       if the TMLE was not run.}
#'     \item{bayes_tmle}{List with posterior draws for \eqn{beta} (\code{samples}), a \code{chains} by \code{draws} by \code{p} array and
#'     the mean acceptance rate (\code{acc_rate}). \code{NULL} if the generalized Bayesian estimator was not run.}
#'   }
#'
#' @section Inference:
#' Standard errors are the empirical standard deviation of the estimated efficient influence function divided by \eqn{\sqrt{n}}.
#' Under the conditions of Theorem 1, the one-step and targeted estimators are asymptotically normal and efficient, with
#' asymptotic variance given by the variance of the efficient influence function.
#'
#' @name automsm
#' @docType class
NULL
