test_that("simulating dose-response example data for works", {
  dat <- simulate_dose_response(N = 100, treatments = 5, sigma = 0.1, seed = 1)
  expect_equal(names(dat), c("X1", "X2", "A", "Y"))
  expect_true(all(dat$A %in% 1:5))
  expect_equal(nrow(dat), 100)
})

test_that("dose-response one-step estimation works on simulated dataset", {
  dat <- simulate_dose_response(N = 500, treatments = 5, sigma = 0.1, seed = 1, nonlinear = FALSE)
  set.seed(1)
  fit <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + A,
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    outcome = "continuous",
    tmle = FALSE
  )

  expect_equal(fit$onestep$est, c(-0.009, 1.002), tolerance = 0.01)
  expect_equal(fit$onestep$se, c(0.011, 0.004), tolerance = 0.01)
})

test_that("dose-response TMLE works on simulated dataset", {
  dat <- simulate_dose_response(N = 500, treatments = 5, sigma = 0.1, seed = 1, nonlinear = FALSE)
  set.seed(1)
  fit <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + A,
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    outcome = "continuous",
    tmle = TRUE
  )

  expect_equal(fit$tmle$est, c(-0.009, 1.002), tolerance = 0.01)
  expect_equal(fit$tmle$se, c(0.011, 0.004), tolerance = 0.01)
})
