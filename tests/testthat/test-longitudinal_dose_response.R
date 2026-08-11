test_that("simulating longitudinal dose-response example data works", {
  dat <- simulate_longitudinal_dose_response(N = 10, tau = 2, seed = 1)
  expect_equal(names(dat), c("L1", "L2", "A1", "A2", "Y"))
  expect_equal(nrow(dat), 10)
})

test_that("longitudinal dose-response one-step estimation works on simulated dataset", {
  tau <- 2
  dat <- simulate_longitudinal_dose_response(N = 500, tau = tau, seed = 1)

  Ls <- lapply(1:tau, \(t) paste0("L", t))
  As <- paste0("A", 1:tau)
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  set.seed(1)
  fit <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    learners_outcome = "SL.glm",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    outer_folds = 5,
    tmle = FALSE
  )

  expect_equal(fit$onestep$est, c(-1.037, 0.591), tolerance = 0.01)
  expect_equal(fit$onestep$se, c(0.263, 0.176), tolerance = 0.01)
})

test_that("longitudinal dose-response TMLE estimation works on simulated dataset", {
  tau <- 2
  dat <- simulate_longitudinal_dose_response(N = 500, tau = tau, seed = 1)

  Ls <- lapply(1:tau, \(t) paste0("L", t))
  As <- paste0("A", 1:tau)
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  set.seed(1)
  fit <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    learners_outcome = "SL.glm",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    outer_folds = 5,
    tmle = TRUE
  )

  expect_equal(fit$tmle$est, c(-0.978, 0.531), tolerance = 0.01)
  expect_equal(fit$tmle$se, c(0.262, 0.175), tolerance = 0.01)
})

test_that("hand-coding working model with summary measures gives exact same results as summary_measures", {
  tau <- 2
  dat <- simulate_longitudinal_dose_response(N = 500, tau = tau, seed = 1)

  Ls <- lapply(1:tau, \(t) paste0("L", t))
  As <- paste0("A", 1:tau)
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  set.seed(1)
  fit <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    learners_outcome = "SL.glm",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = working_model_linear,
    outer_folds = 5,
    tmle = TRUE
  )

  X_regimes <- torch::torch_cat(
  list(torch::torch_ones(c(nrow(regimes), 1)), torch::torch_tensor(rowSums(regimes))$reshape(c(nrow(regimes), 1))),
  dim = 2)

  set.seed(1)
  fit2 <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    learners_outcome = "SL.glm",
    regimes = regimes,
    loss = loss_weighted_sum(loss_cross_entropy_logit),
    working_model = function(beta, X) {
      n <- dim(X)[1]
      preds <- X_regimes$matmul(beta)$reshape(c(1, nrow(regimes))) # 1 x k
      torch::torch_ones(c(n, 1))$matmul(preds) # n x k
    },
    p = 2,
    formula = ~1,
    outer_folds = 5,
    tmle = TRUE
  )

  expect_equal(fit$onestep$est, c(-1.037, 0.591), tolerance = 0.01)
  expect_equal(fit$onestep$se, c(0.263, 0.176), tolerance = 0.01)
  expect_equal(fit$tmle$est, c(-0.978, 0.531), tolerance = 0.01)
  expect_equal(fit$tmle$se, c(0.262, 0.175), tolerance = 0.01)

  expect_equal(fit$onestep$est, fit2$onestep$est, tolerance = 1e-4)
  expect_equal(fit$onestep$se, fit2$onestep$se, tolerance = 1e-4)
  expect_equal(fit$tmle$est, fit2$tmle$est, tolerance = 1e-4)
  expect_equal(fit$tmle$se, fit2$tmle$se, tolerance = 1e-4)
})
