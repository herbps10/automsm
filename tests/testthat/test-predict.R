test_that("predict() is consistent across types and estimators", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    nuisance = nuisance_control(
      learners_trt = "SL.mean",
      learners_outcome = "SL.mean"
    ),
    outcome_type = "continuous",
    tmle = TRUE
  )

  nd <- data.frame(X1 = seq(0, 1, length.out = 7), X2 = seq(1, 0, length.out = 7))
  for(est in c("tmle", "onestep")) {
    p <- predict(fit, nd, estimator = est)
    m <- predict(fit, nd, estimator = est, type = "joint", seed = 1L)

    expect_length(p, nrow(nd))
    expect_equal(dim(m), c(1e3L, nrow(nd)))
    expect_lt(max(abs(colMeans(m) - p)), 0.05 * max(abs(p)) + 1e-3)
    expect_equal(m, predict(fit, nd, estimator = est, type = "joint", seed = 1L), tolerance = 0)

    r <- predict(fit, nd, estimator = est, type = "rvar", seed = 1L)
    expect_equal(length(r), nrow(nd))
    expect_equal(posterior::draws_of(r), m, ignore_attr = TRUE)
  }
})

test_that("bayes predictions come from bayes_tmle$draws, correctly oriented", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- suppressWarnings(cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    nuisance = nuisance_control(
      learners_trt = "SL.mean",
      learners_outcome = "SL.mean"
    ),
    outcome_type = "continuous",
    tmle = TRUE,
    bayes = bayes_control(warmup = 200, draws = 200, seed = 1, chains = 2)
  ))

  nd <- data.frame(X1 = seq(0, 1, length.out = 7), X2 = seq(1, 0, length.out = 7))

  m <- predict(fit, nd, estimator = "bayes", type = "joint")

  bm <- as.matrix(posterior::as_draws_matrix(fit$bayes_tmle$draws))
  expect_equal(nrow(m), nrow(bm))
  expect_equal(m, t(model.matrix(~X2, nd) %*% t(bm)), ignore_attr = TRUE, tolerance = 1e-7)
  expect_equal(predict(fit, nd, estimator = "bayes"), apply(m, 2, median), ignore_attr = TRUE, tolerance = 1e-7)
})
