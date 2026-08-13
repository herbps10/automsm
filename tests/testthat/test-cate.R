# ----- One-step fixtures -----
test_that("cate continous linear one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_continuous_linear.rds"))

  set.seed(1)
  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = FALSE,
    nuisance = fx$nuisance
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-6)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-6)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-6)
})

test_that("cate continous binary one-step fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_binary_linear.rds"))

  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = FALSE,
    nuisance = fx$nuisance
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-6)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-6)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-6)
})

# ----- TMLE fixtures -----
test_that("cate continous linear TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_continuous_linear.rds"))

  set.seed(1)
  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = TRUE,
    nuisance = fx$nuisance
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-6)
  expect_equal(fit$tmle$se, fx$fit_tmle$tmle$se, tolerance = 1e-6)
})

test_that("cate continous binary TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_binary_linear.rds"))

  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = TRUE,
    nuisance = fx$nuisance
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-5)
  expect_equal(fit$tmle$se,  fx$fit_tmle$tmle$se, tolerance = 1e-5)
})

# ----- Bayes TMLE fixtures -----
test_that("cate continous linear Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_continuous_linear.rds"))

  set.seed(1)
  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = TRUE,
    bayes = TRUE,
    bayes_chains = 2,
    bayes_draws = 100,
    nuisance = fx$nuisance
  )

  expect_equal(
    apply(fit$bayes_tmle$samples, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, mean),
    tolerance = 1e-6
  )
  expect_equal(
    apply(fit$bayes_tmle$samples, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, sd),
    tolerance = 1e-6
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-6)
})

test_that("cate continous binary Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_binary_linear.rds"))

  set.seed(1)
  fit <- cate(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "binomial",
    tmle = TRUE,
    bayes = TRUE,
    bayes_chains = 2,
    bayes_draws = 100,
    nuisance = fx$nuisance
  )

  expect_equal(
    apply(fit$bayes_tmle$samples, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, mean),
    tolerance = 1e-5
  )
  expect_equal(
    apply(fit$bayes_tmle$samples, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$samples, 3, sd),
    tolerance = 1e-5
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-6)
})

# ----- Integration tests -----
test_that("cate one-step runs on continuous linear simulated data", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    outcome = "continuous",
    tmle = FALSE
  )

  expect_length(fit$onestep$est, 2)
  expect_true(all(is.finite(fit$onestep$est)))
  expect_true(all(is.finite(fit$onestep$se)))
})

test_that("cate TMLE runs on continuous linear simulated data", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = TRUE
  )


  expect_length(fit$tmle$est, 2)
  expect_true(all(is.finite(fit$tmle$est)))
  expect_true(all(is.finite(fit$tmle$se)))
})

test_that("cate Bayes TMLE runs on continuous linear simulated data", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = TRUE,
    bayes = TRUE,
    bayes_chains = 2,
    bayes_draws = 100
  )

  expect_equal(dim(fit$bayes_tmle$samples), c(2, 100, 2))
  expect_true(all(is.finite(fit$bayes_tmle$samples)))
  expect_true(is.finite(fit$bayes_tmle$acc_rate))
})

