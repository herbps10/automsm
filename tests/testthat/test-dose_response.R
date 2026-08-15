# ----- One-step fixtures -----
test_that("dose_response continous linear one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_continuous_linear.rds"))

  set.seed(1)
  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = FALSE,
    nuisance_estimates = fx$nuisance
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-3)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-3)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-3)
})

test_that("dose_response continous binary one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_binary_linear.rds"))

  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = FALSE,
    nuisance_estimates = fx$nuisance
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-3)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-3)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-3)
})

# ----- TMLE fixtures -----
test_that("dose_response continous linear TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_continuous_linear.rds"))

  set.seed(1)
  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = TRUE,
    nuisance_estimates = fx$nuisance
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se, fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

test_that("dose_response continous binary TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_binary_linear.rds"))

  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = TRUE,
    nuisance_estimates = fx$nuisance
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se,  fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

# ----- Bayes TMLE fixtures -----
test_that("dose_response continous linear Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_continuous_linear.rds"))

  set.seed(1)
  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100
    ),
    nuisance_estimates = fx$nuisance
  )

  expect_equal(
    apply(fit$bayes_tmle$samples, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, mean),
    tolerance = 1e-2
  )
  expect_equal(
    apply(fit$bayes_tmle$samples, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, sd),
    tolerance = 1e-2
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-1)
})

test_that("dose_response continous binary Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_binary_linear.rds"))

  set.seed(1)
  fit <- dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100,
    ),
    nuisance_estimates = fx$nuisance
  )

  expect_equal(
    apply(fit$bayes_tmle$samples, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, mean),
    tolerance = 1e-4
  )
  expect_equal(
    apply(fit$bayes_tmle$samples, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, sd),
    tolerance = 1e-4
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-4)
})

# ----- Integration tests -----
test_that("dose-response one-step runs on continuous linear simulated data", {
  dat <- simulate_dose_response(N = 500, treatments = 5, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + A,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    ),
    outcome = "continuous",
    tmle = FALSE
  )

  expect_length(fit$onestep$est, 2)
  expect_true(all(is.finite(fit$onestep$est)))
  expect_true(all(is.finite(fit$onestep$se)))
})

test_that("dose-response TMLE runs on continuous linear simulated data", {
  dat <- simulate_dose_response(N = 500, treatments = 5, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + A,
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    ),
    outcome = "continuous",
    tmle = TRUE
  )

  expect_length(fit$tmle$est, 2)
  expect_true(all(is.finite(fit$tmle$est)))
  expect_true(all(is.finite(fit$tmle$se)))
})

test_that("dose-response Bayes TMLE runs on continuous linear simulated data", {
  dat <- simulate_dose_response(N = 500, treatments = 5, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + A,
    outcome = "continuous",
    nuisance = nuisance_control(
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm"
    ),
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100,
    )
  )

  expect_equal(dim(fit$bayes_tmle$samples), c(2, 100, 2))
  expect_true(all(is.finite(fit$bayes_tmle$samples)))
  expect_true(is.finite(fit$bayes_tmle$acc_rate))
})
