test_that("simulating cate example data for works", {
  dat <- simulate_cate(N = 10, sigma = 0.1, seed = 1)
  expect_equal(names(dat), c("X1", "X2", "A", "Y"))
  expect_equal(nrow(dat), 10)
})

test_that("cate one-step estimation works on simulated dataset", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X1,
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    outcome = "continuous",
    tmle = FALSE
  )

  expect_equal(fit$onestep$est, c(0.985, 0.01), tolerance = 0.01)
  expect_equal(fit$onestep$se, c(0.0186, 0.036), tolerance = 0.01)
})

test_that("cate TMLE works on simulated dataset", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X1,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = TRUE
  )

  expect_equal(fit$tmle$est, c(1.018, 0.056), tolerance = 0.01)
  expect_equal(fit$tmle$se, c(0.038, 0.069), tolerance = 0.01)
})

test_that("cate Bayesian TMLE works on simulated dataset", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X1,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = TRUE,
    bayes = TRUE,
    bayes_chains = 2,
    bayes_draws = 100
  )

  bayes_mean <- apply(fit$tmle$samples, 3, mean)
  bayes_sd <- apply(fit$tmle$samples, 3, sd)

  expect_equal(bayes_mean, c(1.020, 0.044), tolerance = 0.01)
  expect_equal(bayes_sd, c(0.093, 0.183), tolerance = 0.01)
})

