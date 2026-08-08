test_that("simulating example data for treatment effect modification works", {
  dat <- simulate_treatment_effect_modification(N = 10, sigma = 0.1, seed = 1)
  expect_equal(names(dat), c("X1", "X2", "A", "Y"))
  expect_equal(nrow(dat), 10)
})

test_that("one-step estimation works on simulated dataset", {
  dat <- simulate_treatment_effect_modification(N = 100, sigma = 0.1, seed = 1)
  set.seed(1)
  fit <- treatment_effect_modification(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X1,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = FALSE
  )

  expect_equal(fit$onestep$est, c(0.02155635, 0.02059106), tolerance = 0.001)
  expect_equal(fit$onestep$se, c(0.2153977, 0.3700055), tolerance = 0.001)
})

test_that("TMLE works on simulated dataset", {
  dat <- simulate_treatment_effect_modification(N = 100, sigma = 0.1, seed = 1)
  set.seed(1)
  fit <- treatment_effect_modification(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X1,
    learners_trt = "SL.mean",
    learners_outcome = "SL.mean",
    tmle = TRUE
  )

  expect_equal(fit$tmle$est, c(0.02168154, 0.01771841), tolerance = 0.01)
  expect_equal(fit$tmle$se, c(0.2135827, 0.3699576), tolerance = 0.01)
})
