#' Validate arguments for longitudinal  NP-MSM estimators
#' @noRd
validate_longitudinal_arguments <- function(
    data, Ls, As, Y, formula, summary_measures,
    outcome_type, loss, working_model, p,
    learners_trt, learners_outcome,
    outer_folds, inner_folds,
    tmle, tmle_maxiter, tmle_linear,
    epsilon, nuisance) {
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)
  checkmate::assert_character(As, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(As, choices = names(data))
  checkmate::assert_list(Ls, len = length(As))
  for(Lt in Ls) checkmate::assert_subset(Lt, choices = names(data))
  checkmate::assert_string(Y)
  checkmate::assert_choice(Y, choices = names(data))
  checkmate::assert_formula(formula)
  checkmate::assert_function(summary_measures, null.ok = TRUE)
  checkmate::assert_choice(outcome_type, choices = c("continuous", "binomial"))
  checkmate::assert_function(loss)
  checkmate::assert_function(working_model)
  checkmate::assert_character(learners_trt, min.len = 1, any.missing = FALSE)
  checkmate::assert_character(learners_outcome, min.len = 1, any.missing = FALSE)
  checkmate::assert_count(outer_folds, positive = TRUE)
  checkmate::assert_count(inner_folds, positive = TRUE)
  checkmate::assert_flag(tmle)
  checkmate::assert_count(tmle_maxiter, positive = TRUE)
  checkmate::assert_flag(tmle_linear)
  checkmate::assert_number(epsilon, lower = 0, upper = 0.5, finite = TRUE)
  checkmate::assert(
    checkmate::check_null(nuisance),
    checkmate::check_list(nuisance),
    combine = "or"
  )
}

#' Validate common arguments for single-time-point NP-MSM estimators
#'
#' @noRd
validate_msm_arguments <- function(
    data, X, A, Y, formula,
    outcome_type, loss, working_model, p,
    learners_trt, learners_outcome,
    outer_folds, inner_folds,
    tmle, tmle_maxiter, tmle_linear,
    bayes, bayes_draws, bayes_chains, bayes_prior,
    epsilon, nuisance) {
  # Argument checks
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)
  checkmate::assert_string(A)
  checkmate::assert_subset(A, choices = names(data))
  checkmate::assert_character(X, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(X, choices = names(data))
  checkmate::assert_string(Y)
  checkmate::assert_choice(Y, choices = names(data))
  checkmate::assert_formula(formula)
  checkmate::assert_choice(outcome_type, choices = c("continuous", "binomial"))
  checkmate::assert_function(loss)
  checkmate::assert_function(working_model)
  checkmate::assert_count(p, positive = TRUE, null.ok = TRUE)
  checkmate::assert_character(learners_trt, min.len = 1, any.missing = FALSE)
  checkmate::assert_character(learners_outcome, min.len = 1, any.missing = FALSE)
  checkmate::assert_count(outer_folds, positive = TRUE)
  checkmate::assert_count(inner_folds, positive = TRUE)
  checkmate::assert_flag(tmle)
  checkmate::assert_count(tmle_maxiter, positive = TRUE)
  checkmate::assert_flag(tmle_linear)
  checkmate::assert_number(epsilon, lower = 0, upper = 0.5, finite = TRUE)

  checkmate::assert(
    checkmate::check_null(nuisance),
    checkmate::check_list(nuisance),
    combine = "or"
  )

  TRUE
}


#' Shared plug-in and one-step estimation
#'
#' @noRd
estimate_plugin_and_onestep <- function(
  Lm_fn, psi, design_matrix, Q, Delta, p, joint_draws, seed
) {
  n <- dim(design_matrix)[1]

  plugin <- B(Lm_fn, psi, design_matrix, Q, p)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  onestep_est <- onestep(
    Lm_fn,
    psi,
    beta,
    design_matrix,
    Q,
    Delta,
    p
  )

  if(!is.null(seed)) set.seed(seed)
  onestep_joint <- mvtnorm::rmvnorm(
    joint_draws,
    mean = as.numeric(onestep_est$est),
    sigma = var(as.matrix(onestep_est$eif)) / n
  )

  list(
    plugin = list(est = as.numeric(plugin)),
    onestep = list(
      est = as.numeric(onestep_est$est),
      se = as.numeric(onestep_est$se),
      lower = as.numeric(onestep_est$lower),
      upper = as.numeric(onestep_est$upper),
      eif = onestep_est$eif,
      joint_draws = onestep_joint
    )
  )
}
