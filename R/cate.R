#' Estimate Marginal Structural Model for the Conditional Average Treatment Effect (CATE)
#'
#' @inheritParams automsm_shared_params
#' @param X A character vector giving the covariate column names in \code{data}.
#' @param A A string naming the treatment column in \code{data}.
#' @param Y A string naming the outcome column in \code{data}.
#'
#' @return An object of class \code{"automsm"}; see [automsm] for the full list of components.
#'
#' @seealso [dose_response()] for a categorical point treatment,
#'   [longitudinal_dose_response()] for time-varying treatment.
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
  nuisance = nuisance_control(),
  nuisance_estimates = NULL,
  tmle = tmle_control(),
  bayes = FALSE,
  onestep = onestep_control()
) {
  nuisance <- as_nuisance_control(nuisance)
  tmle <- resolve_fluctuation(as_tmle_control(tmle), outcome_type)
  bayes <- as_bayes_control(bayes)
  onestep <- as_onestep_control(onestep)

  validate_msm_arguments(
    data, X, A, Y, formula, outcome_type, loss, working_model, p, nuisance_estimates
  )

  if (is.null(nuisance_estimates)) {
    nuisance_estimates <- estimate_cate_nuisance(
      data,
      X,
      A,
      Y,
      nuisance,
      outcome_type,
      estimate_conditional_variance = isTRUE(bayes$enabled) && isTRUE(tmle$linear)
    )
  }

  problem <- cate_problem(data, A, Y, formula, p, outcome_type, loss, working_model, nuisance_estimates)

  spec <- msm_spec_cate(tmle_linear = tmle$linear)

  res <- fit_msm(
    problem, spec,
    tmle = tmle,
    bayes = bayes,
    onestep = onestep
  )

  new_automsm(problem, res$base, res$tmle, res$bayes, nuisance)
}

#' @noRd
cate_problem <- function(data, A, Y, formula, p = NULL,
                         outcome_type = "binomial",
                         loss = loss_squared_error,
                         working_model = working_model_linear,
                         nuisance_estimates) {
  n <- nrow(data)
  dm <- build_design_tensor(formula, data, K = 1L, mutate = NULL, p = p)
  pi <- nuisance_estimates$pi
  new_msm_problem(
    estimand = "cate", K = 1L, d = dm$d, p = dm$p, tau = 1L,
    design_matrix = dm$design_matrix,
    Q0 = torch::torch_tensor(rep(1 / n, n)),
    Yt = data[[Y]],
    Lm_fn = Lm(loss, working_model),
    loss = loss, working_model = working_model,
    formula = formula, terms = dm$terms,
    outcome_type = outcome_type,
    nuisance_estimates = nuisance_estimates,
    aux = list(
      A_obs = data[[A]],
      H = as_float_tensor(ifelse(data[[A]] == 1, 1 / pi, -1 / (1 - pi))),
      H0 = as_float_tensor(-1 / (1 - pi)),
      H1 = as_float_tensor(1 / pi)
    )
  )
}

#' @noRd
estimate_cate_nuisance <- function(
  data,
  X,
  A,
  Y,
  control,
  outcome_type,
  estimate_conditional_variance = FALSE
) {
  nuis <- with_nuisance_seed(control, estimate_point_nuisance(
    data, X, A, Y,
    control,
    outcome_type = outcome_type,
    levels = c(0, 1),
    propensity = "binary",
    estimate_conditional_variance = estimate_conditional_variance
  ))

  list(
    pi       = nuis$pi_a[, 2],
    mu0      = nuis$mu_a[, 1],
    mu1      = nuis$mu_a[, 2],
    mu       = nuis$mu,
    condvar0 = nuis$condvar_a[, 1],
    condvar1 = nuis$condvar_a[, 2],
    condvar  = nuis$condvar
  )
}


