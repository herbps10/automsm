#' Estimate Marginal Structural Model for the Conditional Average Treatment Effect (CATE)
#'
#' @param data A \code{data.frame} containing the columns referenced by \code{X},
#'   \code{A}, and \code{Y}.
#' @param X A character vector giving the covariate column names in \code{data}.
#' @param A A string naming the treatment column in \code{data}.
#' @param Y A string naming the outcome column in \code{data}.
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
#' @param learners_trt A character vector of \pkg{SuperLearner} libraries for estimating
#'   the propensity scores.
#' @param learners_outcome A character vector of \pkg{SuperLearner} libraries for
#'   estimating the outcome regressions.
#' @param outer_folds Number of folds in the outer cross-fitting loop.
#' @param inner_folds Number of folds for the inner \pkg{SuperLearner} cross-validation
#'   within each outer cross-fitting fold.
#' @param tmle Logical; whether to run the TMLE estimator.
#' @param tmle_maxiter Maximum number of TMLE iterations.
#' @param tmle_linear Logical; whether to use a linear TMLE fluctuation model
#'   (\code{TRUE}) or a logistic fluctuation model (\code{FALSE}).
#' @param bayes whether to run Bayesian TMLE estimator
#' @param bayes_draws number of MCMC samples
#' @param bayes_chains number of independent MCMC chains
#' @param bayes_prior prior to apply to marginal structural model parameters
#' @param epsilon Adjustment bounding estimated propensities/means away from 0 and 1.
#' @param nuisance Optional list of pre-computed nuisance parameters. If \code{NULL}
#'   (the default), nuisance parameters are estimated internally via cross-fitting.
#'
#' #' @return An object of class \code{"automsm"}: a list with components
#'   \describe{
#'     \item{estimand}{Character string, \code{"cate"}.}
#'     \item{p}{Number of working-model coefficients.}
#'     \item{n}{Sample size.}
#'     \item{formula}{The working-model formula used.}
#'     \item{working_model, loss}{The working model and loss function used.}
#'     \item{terms}{Character vector of working-model design-matrix term names.}
#'     \item{learners_trt, learners_outcome}{The \pkg{SuperLearner} libraries used.}
#'     \item{nuisance}{The (estimated or supplied) nuisance parameters.}
#'     \item{plugin}{A list with the plug-in point estimate (\code{est}).}
#'     \item{onestep}{A list with the one-step point estimate (\code{est}), standard
#'       errors (\code{se}), confidence-interval bounds (\code{lower}, \code{upper}),
#'       the estimated efficient influence function (\code{eif}), and joint draws
#'       (\code{joint_draws}).}
#'   }

#' @importFrom torch torch_tensor torch_zeros torch_reshape nn_bce_with_logits_loss nn_mse_loss
#' @importFrom stats dnorm qnorm
#' @importFrom adaptMCMC MCMC
#'
#' @export
cate <- function(
  data,
  X,
  A,
  Y,
  formula,
  outcome_type = "binomial",
  loss = loss_squared_error,
  working_model = working_model_linear,
  p = NULL,
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm",
  outer_folds = 5,
  inner_folds = 5,
  tmle = TRUE,
  tmle_maxiter = 25,
  tmle_linear = TRUE,
  bayes = FALSE,
  bayes_draws = 1e3,
  bayes_chains = 4,
  bayes_prior = function(beta) {
    sum(stats::dnorm(as.numeric(beta), mean = 0, sd = 1, log = TRUE))
  },
  epsilon = 1e-5,
  nuisance = NULL
) {
  validate_msm_arguments(
    data, X, A, Y, formula,
    outcome_type, loss, working_model, p,
    learners_trt, learners_outcome,
    outer_folds, inner_folds,
    tmle, tmle_maxiter, tmle_linear,
    bayes, bayes_draws, bayes_chains, bayes_prior,
    epsilon, nuisance
  )

  n <- nrow(data)

  if (is.null(nuisance)) {
    nuisance <- estimate_cate_nuisance(
      data,
      X,
      A,
      Y,
      learners_trt,
      learners_outcome,
      outer_folds,
      inner_folds,
      outcome_type,
      estimate_conditional_variance = bayes && tmle_linear == TRUE,
      epsilon = epsilon
    )
  }

  dm <- build_design_tensor(formula, data, K = 1L, mutate = NULL, p = p)
  pi <- nuisance$pi

  problem <- new_msm_problem(
    estimand = "cate", K = 1L, d = dm$d, p = dm$p, tau = 1L,
    design_matrix = dm$design_matrix,
    Q0 = torch::torch_tensor(rep(1 / n, n)),
    Yt = data[[Y]],
    Lm_fn = Lm(loss, working_model),
    loss = loss, working_model = working_model,
    formula = formula, terms = dm$terms,
    outcome_type = outcome_type, nuisance = nuisance,
    aux = list(
      H = as_float_tensor(ifelse(data[[A]] == 1, 1 / pi, -1 / (1 - pi))),
      H0 = as_float_tensor(-1 / (1 - pi)),
      H1 = as_float_tensor(1 / pi)
    )
  )

  spec <- msm_spec_cate(tmle_linear = tmle_linear)

  res <- fit_msm(
    problem, spec,
    tmle = tmle_control(tmle, tmle_maxiter, tmle_linear),
    bayes = bayes_control(bayes, bayes_draws, bayes_chains, bayes_prior),
    onestep = onestep_control(1e3, NULL)
  )

  assemble_result(
    "cate", problem$p, problem$d, problem$n, formula, working_model, loss, problem$terms,
    learners_trt, learners_outcome, nuisance,
    plugin = res$base$plugin, onestep = res$base$onestep,
    tmle = res$tmle, bayes_tmle = res$bayes, tau = 1L
  )
}

#' @noRd
estimate_cate_nuisance <- function(
  data,
  X,
  A,
  Y,
  learners_trt,
  learners_outcome,
  outer_folds,
  inner_folds,
  outcome_type,
  estimate_conditional_variance = FALSE,
  epsilon = 1e-5
) {
  nuis <- estimate_point_nuisance(
    data, X, A, Y,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    outcome_type = outcome_type,
    levels = c(0, 1),
    propensity = "binary",
    estimate_conditional_variance = estimate_conditional_variance,
    epsilon = epsilon
  )

  list(
    pi      = nuis$pi_a[, 2],
    mu0     = nuis$mu_a[, 1],
    mu1     = nuis$mu_a[, 2],
    mu      = nuis$mu,
    condvar = nuis$condvar
  )
}


