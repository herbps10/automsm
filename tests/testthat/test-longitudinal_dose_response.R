# ----- One-step fixtures -----
test_that("dose_response continous one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "longitudinal_dose_response_continuous.rds"))

  set.seed(1)
  fit <- longitudinal_dose_response(
    fx$dat,
    fx$Ls,
    fx$As,
    "Y",
    regimes = fx$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_squared_error,
    outcome_type = "continuous",
    tmle = FALSE,
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-3)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-3)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-3)
})

test_that("dose_response binary one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "longitudinal_dose_response_binary.rds"))

  fit <- longitudinal_dose_response(
    fx$dat,
    Ls = fx$Ls,
    As = fx$As,
    "Y",
    regimes = fx$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_cross_entropy_logit,
    outcome_type = "binomial",
    tmle = FALSE,
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-3)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-3)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-3)
})

# ----- TMLE fixtures -----
test_that("dose_response continous linear TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "longitudinal_dose_response_continuous.rds"))

  fit <- longitudinal_dose_response(
    fx$dat,
    fx$Ls,
    fx$As,
    "Y",
    regimes = fx$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_squared_error,
    outcome_type = "continuous",
    tmle = TRUE,
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se, fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

test_that("dose_response continous binary TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "longitudinal_dose_response_binary.rds"))

  fit <- longitudinal_dose_response(
    fx$dat,
    fx$Ls,
    fx$As,
    "Y",
    regimes = fx$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_cross_entropy_logit,
    outcome_type = "binomial",
    tmle = TRUE,
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se,  fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

# ----- Integration tests -----

integration_test_setup <- function(N = 500, tau = 2) {
  dat <- simulate_longitudinal_dose_response(N = N, tau = tau, seed = 1)

  Ls <- lapply(1:tau, \(t) paste0("L", t))
  As <- paste0("A", 1:tau)
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  list(data = dat, Ls = Ls, As = As, regimes = regimes, tau = tau)
}

test_that("longitudinal dose-response one-step estimation runs on continuous simulated data", {
  setup <- integration_test_setup()

  set.seed(1)
  fit <- longitudinal_dose_response(
    setup$data,
    setup$Ls,
    setup$As,
    "Y",
    regimes = setup$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    tmle = FALSE,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    )
  )

  expect_length(fit$onestep$est, 2)
  expect_true(all(is.finite(fit$onestep$est)))
  expect_true(all(is.finite(fit$onestep$se)))
})

test_that("longitudinal dose-response TMLE estimation runs on continuous simulated data", {
  setup <- integration_test_setup()

  set.seed(1)
  fit <- longitudinal_dose_response(
    setup$data,
    setup$Ls,
    setup$As,
    "Y",
    regimes = setup$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    tmle = TRUE,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    )
  )

  expect_length(fit$tmle$est, 2)
  expect_true(all(is.finite(fit$tmle$est)))
  expect_true(all(is.finite(fit$tmle$se)))
})

test_that("hand-coding working model with summary measures gives exact same results as summary_measures", {
  setup <- integration_test_setup()

  set.seed(1)
  fit <- longitudinal_dose_response(
    setup$data,
    setup$Ls,
    setup$As,
    "Y",
    regimes = setup$regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    tmle = TRUE,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    )
  )

  X_regimes <- torch::torch_cat(
    list(
      torch::torch_ones(c(nrow(setup$regimes), 1)),
      torch::torch_tensor(rowSums(setup$regimes))$reshape(c(nrow(setup$regimes), 1))
    ),
    dim = 2
  )

  set.seed(1)
  fit2 <- longitudinal_dose_response(
    setup$data,
    setup$Ls,
    setup$As,
    "Y",
    regimes = setup$regimes,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = function(beta, X) {
      n <- dim(X)[1]
      preds <- X_regimes$matmul(beta)$reshape(c(1, nrow(setup$regimes))) # 1 x k
      torch::torch_ones(c(n, 1))$matmul(preds) # n x k
    },
    p = 2,
    formula = ~1,
    tmle = TRUE,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    )
  )

  expect_equal(fit$onestep$est, fit2$onestep$est, tolerance = 1e-4)
  expect_equal(fit$onestep$se, fit2$onestep$se, tolerance = 1e-4)
  expect_equal(fit$tmle$est, fit2$tmle$est, tolerance = 1e-4)
  expect_equal(fit$tmle$se, fit2$tmle$se, tolerance = 1e-4)
})

run_longitudinal_regression <- function(setup) {
  control <- nuisance_control()
  cv <- nuisance_setup(nrow(setup$data), control, outcome_type = "binomial")

  estimate_longitudinal_dose_response_regressions(
    data = setup$data, Ls = setup$Ls, As = setup$As, Y = "Y",
    regimes = setup$regimes, learners_outcome = control$learners_outcome,
    cv = cv$cv, outcome_type = "binomial", outcome_family = cv$outcome_family, cv_control = cv$cv_control,
    epsilon = control$epsilon
  )
}

test_that("regression output has shape k x n x (tau + 1) and is fully populated", {
  set.seed(1)
  setup <- integration_test_setup()
  res <- run_longitudinal_regression(setup)

  expect_true(is.array(res))
  expect_equal(dim(res), c(nrow(setup$regimes), nrow(setup$data), setup$tau + 1L))
  expect_false(anyNA(res))
  expect_true(all(is.finite(res)))
})

test_that("regression output has node tau + 1 as the observed outcome, identically for every regime", {
  set.seed(1)
  setup <- integration_test_setup()
  res <- run_longitudinal_regression(setup)

  for(j in seq_len(nrow(setup$regimes))) {
    expect_equal(res[j, , setup$tau + 1L], as.numeric(setup$data$Y))
  }
})

test_that("regression output has slices aligned to the rows of regimes", {
  setup <- integration_test_setup()
  set.seed(1)
  res1 <- run_longitudinal_regression(setup)

  perm <- c(3L, 1L, 4L, 2L)
  regimes_perm <- setup$regimes[perm, , drop = FALSE]
  rownames(setup$regimes) <- NULL
  setup$regimes <- regimes_perm

  set.seed(1)
  res2 <- run_longitudinal_regression(setup)

  for(j in seq_along(perm)) {
    expect_equal(res2[j, , ], res1[perm[j], , ], tolerance = 1e-7)
  }
})

test_that("regression output are identical for duplicate regimes", {
  setup <- integration_test_setup()
  setup$regimes <- rbind(setup$regimes, setup$regimes)
  res <- run_longitudinal_regression(setup)

  expect_equal(dim(res)[1], nrow(setup$regimes))
  k <- nrow(setup$regimes) / 2
  for(j in seq_len(k)) {
    expect_equal(res[j,,], res[k + j, ,], tolerance = 1e-7)
  }
})

