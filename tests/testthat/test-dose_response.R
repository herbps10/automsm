source(test_path("test-internals.R"))
source(test_path("helper-invariants.R"))
source(test_path("helper-bayes.R"))

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
    nuisance_estimates = fx$nuisance_estimates
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
    nuisance_estimates = fx$nuisance_estimates
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
    nuisance_estimates = fx$nuisance_estimates
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
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se,  fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

# ----- Bayes TMLE fixtures -----
test_that("dose_response continous linear Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_continuous_linear.rds"))

  set.seed(1)
  fit <- suppressWarnings(dose_response(
    fx$dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = "continuous",
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100,
      warmup = 100
    ),
    nuisance_estimates = fx$nuisance_estimates
  ))

  expect_equal(
    apply(fit$bayes_tmle$draws, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, mean),
    tolerance = 1e-2
  )
  expect_equal(
    apply(fit$bayes_tmle$draws, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, sd),
    tolerance = 1e-2
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-1)
})

test_that("dose_response continous binary Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "dose_response_binary_linear.rds"))

  set.seed(1)
  fit <- suppressWarnings(dose_response(
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
      warmup = 100,
      fluctuation = "logistic"
    ),
    nuisance_estimates = fx$nuisance
  ))

  expect_equal(
    apply(fit$bayes_tmle$draws, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, mean),
    tolerance = 1e-4
  )
  expect_equal(
    apply(fit$bayes_tmle$draws, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, sd),
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
  dat <- simulate_dose_response(N = 500, treatments = 6, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- suppressWarnings(dose_response(
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
      chains = 1,
      draws = 100,
      warmup = 200
    )
  ))

  expect_equal(dim(fit$bayes_tmle$draws), c(100, 1, 2))
  expect_true(all(is.finite(fit$bayes_tmle$draws)))
  expect_true(is.finite(fit$bayes_tmle$acc_rate))
})

test_that("dose-response Bayes TMLE BvM holds with multiple treatments", {
  dat <- simulate_dose_response(N = 500, treatments = 6, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- suppressWarnings(dose_response(
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
      draws = 500,
      warmup = 500
    )
  ))

  s <- fit$bayes_tmle$diagnostics

  expect_lt(max(abs(s$median - fit$tmle$est) / fit$tmle$se), 0.5)
  expect_lt(max(abs(log(s$sd / fit$tmle$se))), log(1.5))
})

test_that("dose_response gradient of the log-loss equals the summed EIF", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_dose_response(data)
  for(tmle_linear in c(FALSE, TRUE)) {
    it <- dose_response_internals(data, nu, formula = ~A, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
    logf <- bayes_loglik_at(it$problem, it$spec, it$fit, it$condvar, tmle_linear)
    g <- autograd_gradient(logf, torch::torch_tensor(rep(0, it$problem$p)))
    se <- it$tmle$se * it$problem$n
    expect_lt(max(abs(abs(g) - abs(colSums(it$tmle$eif))) / se), 1e-3)
  }
})

test_that("dose_response eif rows are reproducible from the clever covariates", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_dose_response(data)
  for(tmle_linear in c(FALSE, TRUE)) {
    it <- dose_response_internals(data, nu, formula = ~A, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
    st <- it$fit$state
    fn <- it$fit$final
    cl <- as.matrix(fn$clever$obs)
    K <- as.matrix(fn$K_Q)
    resid <- as.numeric(it$problem$Yt) - as.numeric(st$mu$obs)
    expect_equal(abs(as.matrix(it$tmle$eif)), abs(resid * cl + K), tolerance = 1e-4)
  }
})

test_that("(B2) holds exactly against the conditional variance for dose_response", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_dose_response(data)
  tmle_linear <- FALSE
  it <- dose_response_internals(data, nu, formula = ~1 + A, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
  logf <- bayes_loglik_at(it$problem, it$spec, it$fit, it$condvar, linear = tmle_linear)

  H <- hessian_from_gradient(logf, rep(0, it$problem$p))
  expect_lt(attr(H, "asymmetry"), 1e-4)

  expect_equal(-H, Lhat_theory(it), tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("(B2) holds statistically against the empirical EIF for dose_response", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_dose_response(data)
  it <- dose_response_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = "logistic"))
  n <- it$problem$n
  A <- Lhat_theory(it)
  B <- crossprod(as.matrix(it$tmle$eif))
  R <- solve(chol(A))
  R <- t(R) %*% B %*% R
  expect_lt(max(abs(R - diag(nrow(R)))), 10 / sqrt(n))
})

test_that("dose_response generalized posterior with binary outcome and linear or logistic fluctuation concentrates as predicted by BvM", {
  data <- sim1_data(n = 1e3)
  for(bayes_fluctuation in c("linear", "logistic")) {
    suppressWarnings(fit <- dose_response(
      data,
      X = paste0("X", 1:4), A = "A", Y = "Y",
      formula = ~1 + A,
      nuisance_estimates = oracle_nuisance_dose_response(data),
      outcome_type = "binomial",
      tmle = TRUE,
      bayes = bayes_control(
        chains = 2,
        warmup = 500,
        draws = 500,
        seed = 101,
        fluctuation = bayes_fluctuation
      )
    ))

    s <- fit$bayes_tmle$diagnostics

    expect_lt(max(abs(s$median - fit$tmle$est) / fit$tmle$se), 0.5)
    expect_lt(max(abs(log(s$sd / fit$tmle$se))), log(1.5))
  }
})

test_that("dose_response generalized posterior with continuous outcome concentrates as predicted by BvM", {
  data <- sim2_data(n = 500L)
  suppressWarnings(fit <- dose_response(
    data,
    X = paste0("X", 1:4), A = "A", Y = "Y",
    formula = ~A,
    nuisance_estimates = oracle_nuisance_dose_response(data),
    outcome_type = "continuous",
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      warmup = 500,
      draws = 1e3,
      seed = 101
    )
  ))

  s <- fit$bayes_tmle$diagnostics

  expect_lt(max(abs(s$median - fit$tmle$est) / fit$tmle$se), 0.5)
  expect_lt(max(abs(log(s$sd / fit$tmle$se))), log(1.5))
})
