# These helpers return the problem/spec/plugin-in state/fit rather than the
# assembled automsm object.

msm_internals <- function(problem, spec, tmle, bayes, onestep, joint_seed = 1L) {
  res <- automsm:::fit_msm(problem, spec, tmle = tmle, bayes = bayes, onestep = onestep)

  # Precompute psi / Minv / K_Q / clever for use in tests
  state0 <- res$state0
  psi0 <- spec$psi_from_state(problem, state0)
  Minv0 <- automsm:::normalizing_matrix(problem$Lm_fn, psi0, state0$beta, problem$design_matrix, state0$Q, problem$p)

  structure(list(
    problem = problem,
    spec = spec,
    control = list(tmle = tmle, bayes = bayes, onestep = onestep),
    base = res$base,
    plugin = res$base$plugin$est,
    state0 = state0,
    psi0 = psi0,
    Minv0 = Minv0,
    K_Q0 = automsm:::calculate_K(problem$Lm_fn, psi0, state0$beta, problem$design_matrix, Minv0, problem$batched_beta),
    clever0 = spec$make_clever(problem, spec$steps(problem)[[1]], state0, psi0, Minv0),
    fit = res$fit,
    tmle = res$tmle,
    bayes = res$bayes,
    prior = bayes$prior,
    condvar = if(is.null(problem$nuisance$condvar)) NULL else automsm:::as_float_tensor(problem$nuisance$condvar),
    tmle_linear = tmle$linear
  ), class = "msm_internals")
}

cate_internals <- function(
    data, nuisance_estimates, formula = ~X4, A = "A", Y = "Y", p = NULL,
    outcome_type = "binomial", loss = loss_squared_error, working_model = working_model_linear, tmle = tmle_control(), bayes = FALSE, onestep = onestep_control(), ...
) {
  tmle <- automsm:::resolve_fluctuation(automsm:::as_tmle_control(tmle), outcome_type)
  bayes <- automsm:::as_bayes_control(bayes)
  onestep <- automsm:::as_onestep_control(onestep)
  problem <- automsm:::cate_problem(
    data, A = A, Y = Y, formula = formula, p = p, outcome_type = outcome_type, loss = loss,
    working_model = working_model, nuisance_estimates = nuisance_estimates
  )

  msm_internals(problem, automsm:::msm_spec_cate(tmle_linear = tmle$linear), tmle, bayes, onestep, ...)
}

dose_response_internals <- function(
    data, nuisance_estimates, formula = ~A, A = "A", Y = "Y", p = NULL,
    outcome_type = "binomial", loss = loss_squared_error, working_model = working_model_linear, tmle = tmle_control(), bayes = FALSE, onestep = onestep_control(), ...
) {
  tmle <- automsm:::resolve_fluctuation(automsm:::as_tmle_control(tmle), outcome_type)
  bayes <- automsm:::as_bayes_control(bayes)
  onestep <- automsm:::as_onestep_control(onestep)
  problem <- automsm:::dose_response_problem(
    data, A = A, Y = Y, formula = formula, p = p, outcome_type = outcome_type, loss = loss,
    working_model = working_model, nuisance_estimates = nuisance_estimates
  )

  msm_internals(problem, automsm:::msm_spec_dose_response(tmle_linear = tmle$linear), tmle, bayes, onestep, ...)
}

source(test_path("helper-invariants.R"))

test_that("internals helpers agree with exported wrappers", {
  data <- sim1_data(n = 150L)
  nu <- oracle_nuisance_cate(data)

  exported <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y", formula = ~X4,
                   tmle = tmle_control(), bayes = FALSE, nuisance_estimates = nu)
  it <- cate_internals(data, nu, formula = ~X4, tmle = TRUE, joint_seed = 1L)

  expect_equal(it$plugin, exported$plugin$est, tolerance = 0)
  expect_equal(it$base$onestep$est, exported$onestep$est, tolerance = 0)
  expect_equal(it$tmle$est, exported$tmle$est, tolerance = 0)
  expect_equal(it$problem$p, exported$p)
  expect_equal(it$problem$d, exported$d)
  expect_equal(it$problem$terms, exported$terms)
})

