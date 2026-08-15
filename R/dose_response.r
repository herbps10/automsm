#' Non-parametric marginal structural model estimation for dose response curves.
#'
#' @inheritParams automsm_shared_params
#' @param X A character vector giving the covariate column names in \code{data}.
#' @param A A string naming the treatment column in \code{data}.
#' @param Y A string naming the outcome column in \code{data}.
#'
#' @return An object of class \code{"automsm"}; see [automsm] for the full list of components.
#'
#' @seealso [cate()] for conditional average treatment effects,
#'   [longitudinal_dose_response()] for time-varying treatment.
#'
#' @export
dose_response <- function(
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
  tmle <- as_tmle_control(tmle)
  bayes <- as_bayes_control(bayes)
  onestep <- as_onestep_control(onestep)

  validate_msm_arguments(
    data, X, A, Y,
    formula, outcome_type, loss, working_model, p,
    nuisance_estimates
  )

  if (is.null(nuisance_estimates)) {
    nuisance_estimates <- estimate_dose_response_nuisance(
      data,
      X,
      A,
      Y,
      nuisance,
      outcome_type,
      estimate_conditional_variance = isTRUE(bayes$enabled) && isTRUE(tmle$linear)
    )
  }

  problem <- dose_response_problem(data, A, Y, formula, p, outcome_type, loss, working_model, nuisance_estimates)

  spec <- msm_spec_dose_response(tmle_linear = tmle$linear)

  res <- fit_msm(
    problem, spec,
    tmle = tmle,
    bayes = bayes,
    onestep = onestep
  )

  new_automsm(problem, res$base, res$tmle, res$bayes, nuisance)
}

#' @noRd
dose_response_problem <- function(data, A, Y, formula, p = NULL,
                                  outcome_type = "binomial",
                                  loss = loss_squared_error,
                                  working_model = working_model_linear,
                                  nuisance_estimates) {
  n <- nrow(data)
  mat <- stats::model.matrix(formula, data = data)
  terms <- colnames(mat)
  d <- ncol(mat)
  if(is.null(p)) p <- d

  As <- sort(unique(data[[A]]))
  K <- length(As)

  dm <- build_design_tensor(
    formula, data, K = K,
    mutate = function(d, k) {
      d[[A]] <- As[k]
      d
    },
    p = p
  )

  H <- torch_zeros(c(n, K))
  HA <- torch_zeros(c(n, K))
  for (k in 1:K) {
    H[, k] <- (data[[A]] == As[k]) / nuisance_estimates$pi
    HA[, k] <- 1 / nuisance_estimates$pi_a[, k]
  }

  new_msm_problem(
    estimand = "dose_response", K = K, d = dm$p, p = dm$p, tau = 1L,
    design_matrix = dm$design_matrix,
    Q0 = torch::torch_tensor(rep(1 / n, n)),
    Yt = data[[Y]],
    Lm_fn = Lm(loss, working_model),
    loss = loss, working_model = working_model,
    formula = formula, terms = dm$terms,
    outcome_type = outcome_type,
    nuisance_estimates = nuisance_estimates,
    aux = list(
      A_index = match(data[[A]], As),
      H = H,
      HA = HA,
      As = As
    )
  )
}

#' @noRd
estimate_dose_response_nuisance <- function(
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
    propensity = "one_vs_rest",
    estimate_conditional_variance = estimate_conditional_variance
  ))

  list(
    pi_a    = nuis$pi_a,
    mu_a    = nuis$mu_a,
    mu      = nuis$mu,
    pi      = nuis$pi_obs,
    condvar = nuis$condvar
  )
}

