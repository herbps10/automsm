source(test_path("test-internals.R"))
source(test_path("helper-invariants.R"))
source(test_path("helper-bayes.R"))

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
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-6)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-6)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-5)
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
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$plugin$est, fx$fit_onestep$plugin$est, tolerance = 1e-5)
  expect_equal(fit$onestep$est, fx$fit_onestep$onestep$est, tolerance = 1e-5)
  expect_equal(fit$onestep$se, fx$fit_onestep$onestep$se, tolerance = 1e-5)
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
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-5)
  expect_equal(fit$tmle$se, fx$fit_tmle$tmle$se, tolerance = 1e-5)
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
    nuisance_estimates = fx$nuisance_estimates
  )

  expect_equal(fit$tmle$est, fx$fit_tmle$tmle$est, tolerance = 1e-3)
  expect_equal(fit$tmle$se,  fx$fit_tmle$tmle$se, tolerance = 1e-3)
})

# ----- Bayes TMLE fixtures -----
test_that("cate continous linear Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_continuous_linear.rds"))

  set.seed(1)
  fit <- suppressWarnings(cate(
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
      warmup = 100,
      fluctuation = "linear"
    ),
    nuisance_estimates = fx$nuisance_estimates
  ))

  expect_equal(fit$bayes_tmle$fluctuation, "linear")

  expect_equal(
    apply(fit$bayes_tmle$draws, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, mean),
    tolerance = 1e-3
  )
  expect_equal(
    apply(fit$bayes_tmle$draws, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, sd),
    tolerance = 1e-3
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-4)
})

test_that("cate continous binary Bayes TMLE fixture is stable", {
  fx <- readRDS(test_path("fixtures", "cate_binary_linear.rds"))

  set.seed(1)
  fit <- suppressWarnings(cate(
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
    nuisance_estimates = fx$nuisance_estimates
  ))

  expect_equal(fit$bayes_tmle$fluctuation, "logistic")

  expect_equal(
    apply(fit$bayes_tmle$draws, 3, mean),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, mean),
    tolerance = 1e-3
  )
  expect_equal(
    apply(fit$bayes_tmle$draws, 3, sd),
    apply(fx$fit_bayes$bayes_tmle$draws, 3, sd),
    tolerance = 1e-3
  )
  expect_equal(fit$bayes_tmle$acc_rate, fx$fit_bayes$bayes_tmle$acc_rate, tolerance = 1e-4)
})

# ---- Bayes -----
test_that("cate bayes posterior is a well-formed posterior::draws object", {
  data <- sim1_data(n = 200L)
  suppressWarnings(fit <- cate(
    data,
    X = paste0("X", 1:4), A = "A", Y = "Y",
    formula = ~X4,
    nuisance_estimates = oracle_nuisance_cate(data),
    bayes = bayes_control(
      chains = 2,
      draws = 50,
      warmup = 50
    )
  ))

  d <- posterior::as_draws_array(fit)

  # Check the layout
  expect_equal(posterior::niterations(d), 50)
  expect_equal(posterior::nchains(d), 2)
  expect_equal(posterior::nvariables(d), fit$p)
  expect_equal(posterior::variables(d), sprintf("beta[%d]", seq_len(fit$p)))

  expect_s3_class(fit$bayes_tmle$diagnostics, "draws_summary")
  expect_equal(nrow(fit$bayes_tmle$diagnostics), fit$p)
  expect_true(all(is.finite(array(d))))

  rv <- posterior::as_draws_rvars(fit)
  expect_named(rv, "beta")
  expect_equal(length(rv$beta), fit$p)

  expect_equal(nrow(posterior::summarize_draws(fit)), fit$p)

  e <- posterior::as_draws_array(fit, parameter = "epsilon")
  expect_true("lp__" %in% posterior::variables(e))
})

test_that("cate as_draws() errors informatively without a Bayesian fit", {
  data <- sim1_data(n = 100L)

  suppressWarnings(fit <- cate(
    data,
    X = paste0("X", 1:4), A = "A", Y = "Y",
    formula = ~X4,
    nuisance_estimates = oracle_nuisance_cate(data),
    bayes = FALSE
  ))
  expect_error(posterior::as_draws_array(fit), "no generalized Bayesian posterior")
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
    nuisance = nuisance_control(
      learners_trt = "SL.mean",
      learners_outcome = "SL.mean"
    ),
    outcome_type = "continuous",
    tmle = FALSE
  )

  expect_length(fit$onestep$est, 2)
  expect_true(all(is.finite(fit$onestep$est)))
  expect_true(all(is.finite(fit$onestep$se)))
})

test_that("cate TMLE runs on continuous linear simulated data", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
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
    tmle = TRUE
  )

  expect_length(fit$tmle$est, 2)
  expect_true(all(is.finite(fit$tmle$est)))
  expect_true(all(is.finite(fit$tmle$se)))
})

test_that("cate Bayes TMLE runs on continuous linear simulated data", {
  dat <- simulate_cate(N = 500, sigma = 0.1, seed = 1, nonlinear = FALSE, binary = FALSE)
  set.seed(1)
  fit <- suppressWarnings(cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    nuisance = nuisance_control(
      learners_trt = "SL.mean",
      learners_outcome = "SL.mean"
    ),
    formula = ~1 + X2,
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 50,
      warmup = 50
    )
  ))

  expect_equal(dim(fit$bayes_tmle$draws), c(50, 2, 2))
  expect_true(all(is.finite(fit$bayes_tmle$draws)))
  expect_true(is.finite(fit$bayes_tmle$acc_rate))
})

test_that("gradient of the log-loss equals the summed EIF", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  for(tmle_linear in c(FALSE, TRUE)) {
    it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
    logf <- bayes_loglik_at(it$problem, it$spec, it$fit, it$condvar, linear = tmle_linear)
    g <- autograd_gradient(logf, torch::torch_tensor(rep(0, it$problem$p)))
    se <- it$tmle$se * it$problem$n
    expect_lt(max(abs(abs(g) - abs(colSums(it$tmle$eif))) / se), 1e-3)
  }
})

test_that("eif rows are reproducible from the clever covariates", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  for(tmle_linear in c(FALSE, TRUE)) {
    it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
    st <- it$fit$state
    fn <- it$fit$final
    cl <- as.matrix(fn$clever$obs)
    K <- as.matrix(fn$K_Q)
    resid <- as.numeric(it$problem$Yt) - as.numeric(st$mu$obs)
    expect_equal(abs(as.matrix(it$tmle$eif)), abs(resid * cl + K), tolerance = 1e-4)
  }
})

test_that("(B2) holds exactly against the conditional variance for cate", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  tmle_linear <- FALSE
  it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = if(tmle_linear) "linear" else "logistic"))
  logf <- bayes_loglik_at(it$problem, it$spec, it$fit, it$condvar, linear = tmle_linear)

  H <- hessian_from_gradient(logf, rep(0, it$problem$p))
  expect_lt(attr(H, "asymmetry"), 1e-4)

  expect_equal(-H, Lhat_theory(it), tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("(B2) holds statistically against the empirical EIF for cate", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = "logistic"))
  n <- it$problem$n
  A <- Lhat_theory(it)
  B <- crossprod(as.matrix(it$tmle$eif))
  R <- solve(chol(A))
  R <- t(R) %*% B %*% R
  expect_lt(max(abs(R - diag(nrow(R)))), 10 / sqrt(n))
})

test_that("cate generalized posterior with binary outcome and logistic or linear fluctuations concentrates as predicted by BvM", {
  data <- sim1_data(n = 1e3L)
  for(bayes_fluctuation in c("linear", "logistic")) {
    suppressWarnings(fit <- cate(
      data,
      X = paste0("X", 1:4), A = "A", Y = "Y",
      formula = ~X4,
      nuisance_estimates = oracle_nuisance_cate(data),
      outcome_type = "binomial",
      tmle = tmle_control(fluctuation = "logistic"),
      bayes = bayes_control(
        chains = 2,
        warmup = 500,
        draws = 500,
        seed = 101,
        fluctuation = bayes_fluctuation
      )
    ))

    expect_equal(fit$bayes_tmle$linear, bayes_fluctuation == "linear")
    expect_equal(fit$bayes_tmle$fluctuation, bayes_fluctuation)

    s <- fit$bayes_tmle$diagnostics
    expect_lt(max(abs(s$median - fit$tmle$est) / fit$tmle$se), 0.5)
    expect_lt(max(abs(log(s$sd / fit$tmle$se))), log(1.5))
    expect_true(all(s$rhat < 1.05))
  }
})

test_that("cate generalized posterior with continuous outcome concentrates as predicted by BvM", {
  data <- sim2_data(n = 500L)
  suppressWarnings(fit <- cate(
    data,
    X = paste0("X", 1:4), A = "A", Y = "Y",
    formula = ~X4,
    nuisance_estimates = oracle_nuisance_cate(data),
    outcome_type = "continuous",
    tmle = tmle_control(),
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
  expect_true(all(s$rhat < 1.05))
})

test_that("cate posterior is unchanged by fast path", {
  data <- sim1_data(n = 1e3L)
  for(bayes_fluctuation in c("linear", "logistic")) {
    fit <- function() {
      suppressWarnings(cate(
        data,
        X = paste0("X", 1:4), A = "A", Y = "Y",
        formula = ~X4,
        nuisance_estimates = oracle_nuisance_cate(data),
        outcome_type = "binomial",
        tmle = tmle_control(fluctuation = "logistic"),
        bayes = bayes_control(
          chains = 2,
          warmup = 200,
          draws = 200,
          seed = 101,
          fluctuation = bayes_fluctuation
        )
      ))
    }
    a <- fit()
    b <- withr::with_options(list(automsm.B_wls = FALSE), fit())

    expect_equal(a$bayes_tmle$linear, bayes_fluctuation == "linear")
    expect_equal(a$bayes_tmle$fluctuation, bayes_fluctuation)
    expect_equal(b$bayes_tmle$linear, bayes_fluctuation == "linear")
    expect_equal(b$bayes_tmle$fluctuation, bayes_fluctuation)

    sa <- a$bayes_tmle$diagnostics
    sb <- b$bayes_tmle$diagnostics
    expect_equal(sa$mean, sb$mean, tolerance = 1e-1)
    expect_equal(sa$sd, sb$sd, tolerance = 1e-1)
  }
})

# ------ Newton solver -----
test_that("analytic GLM agrees with spec$mu_loss (value and gradient)", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  for(linear in c(TRUE, FALSE)) {
    it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = if(linear) "linear" else "logistic"))

    problem <- it$problem
    spec <- it$spec
    state <- it$state0

    step <- spec$steps(problem)[[1]]
    psi <- spec$psi_from_state(problem, state)
    Minv <- normalizing_matrix(problem$Lm_fn, psi, state$beta, problem$design_matrix, state$Q, problem$p)
    cl <- spec$make_clever(problem, step, state, psi, Minv)
    cl$offset <- spec$fluctuation_offset(problem, step, state, linear = linear)
    g <- spec$fluctuation_glm(problem, step, state, cl)
    fam <- fluctuation_family(linear)

    for(e in list(rep(0, problem$p), c(0.03, -0.02), c(-0.11, 0.07))) {
      eta <- as.vector(g$offset + g$X %*% e)
      expect_equal(sum(fam$dev(eta, g$target)), as.numeric(spec$mu_loss(torch::torch_tensor(e), problem, step, state, cl, linear = linear)), tolerance = 1e-4)
      et <- torch::torch_tensor(e, requires_grad = TRUE)
      ag <- as.numeric(torch::autograd_grad(spec$mu_loss(et, problem, step, state, cl, linear = linear), et)[[1]])
      expect_equal(as.vector(crossprod(g$X, fam$grad(eta, g$target))), ag, tolerance = 1e-4)
    }
  }
})
