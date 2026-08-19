# Implementation-independent tests
#
# These tests assert mathematical properties that must hold for
# ANY correct implementation.

source(test_path("helper-invariants.R"))
source(test_path("helper-bayes.R"))
source(test_path("test-internals.R"))


# ----- Kernel invariants -----

test_that("Lm composes loss(prediction, target), not loss(target, prediction)", {
  # Lm(t, beta, X) must equal L(m_beta(X), t). This orientation matters
  # for every hand-derived quantity below. nn_mse_loss is symmetric but
  # cross-entropy losses are not.
  n <- 5L; K <- 2L; p <- 3L
  set.seed(1)
  Xa <- array(rnorm(n * K * p), dim = c(n, K, p))
  psi_a <- matrix(rnorm(n * K), n, K)
  b <- rnorm(p)

  Lm_fn <- Lm(loss_squared_error, working_model_linear)
  got <- Lm_fn(torch::torch_tensor(psi_a), leaf(b), torch::torch_tensor(Xa))

  pred <- apply(Xa, c(1, 2), function(x) sum(x * b))
  want <- rowSums((pred - psi_a)^2)
  expect_tensor_close(got, torch::torch_tensor(want))
})

test_that("grad_dL matches the closed form -2 x for linear + squared error", {
  K <- 2L; p <- 3L
  set.seed(2)

  Xa <- array(rnorm(1 * K * p), dim = c(1L, K, p))
  psi_a <- matrix(rnorm(1 * K), 1L, K)
  b <- rnorm(p)

  got <- grad_dL(
    Lm(loss_squared_error, working_model_linear),
    leaf(psi_a), leaf(b), torch::torch_tensor(Xa)
  )

  want <- -2 * t(matrix(Xa[1, , ], nrow = K, ncol = p))
  expect_tensor_close(got, torch::torch_tensor(want))
})

test_that("ddL matches the closed form 2 sum_k x_k x_k' at unit weight", {
  n <- 7L; K <- 2L; p <- 3L
  set.seed(3)
  Xa <- array(rnorm(n * K * p), dim = c(n, K, p))
  psi_a <- matrix(rnorm(n * K), n, K)
  b <- rnorm(p)

  got <- ddL(
    Lm(loss_squared_error, working_model_linear),
    leaf(psi_a), leaf(b), torch::torch_tensor(Xa), p
  )

  want <- matrix(0, p, p)

  for(i in seq_len(n)) {
    for(k in seq_len(K)) {
      x <- Xa[i, k, ]
      want <- want + 2 * outer(x, x)
    }
  }

  expect_tensor_close(got, torch::torch_tensor(want))
})

test_that("normalizing matrix inverts E_Q[Lddot] for any Q", {
  n <- 9L; K <- 2L; p <- 3L
  set.seed(4)

  Xa <- array(rnorm(n * K * p), dim = c(n, K, p))
  psi_a <- matrix(rnorm(n * K), n, K)
  b <- rnorm(p)

  Lm_fn <- Lm(loss_squared_error, working_model_linear)

  check <- function(Qv) {
    Minv <- normalizing_matrix(Lm_fn, leaf(psi_a), leaf(b),
                               torch::torch_tensor(Xa), leaf(Qv), p)
    M2 <- dobjective_dbeta(Lm_fn, leaf(psi_a), leaf(b),
                           torch::torch_tensor(Xa), leaf(Qv), p)

    expect_equal(solve(as.matrix(Minv)), as.matrix(M2$detach()), tolerance = 1e-5)
  }

  check(rep(1/n, n))

  w <- runif(n, 0.5, 2)
  check(w / sum(w))
})

test_that("Q_fluctuation satisfies (M1a) and stays normalized", {
  n <- 8L; p <- 2L
  set.seed(5)
  Kmat <- torch::torch_tensor(matrix(rnorm(n * p), n, p))
  eps0 <- torch::torch_tensor(rep(0, p))

  Qu <- torch::torch_tensor(rep(1/n, n))
  expect_equal(as.numeric(Q_fluctuation(eps0, Kmat, Qu)$sum()), 1, tolerance = 1e-6)

  expect_tensor_close(Q_fluctuation(eps0, Kmat, Qu), torch::torch_tensor(rep(1/n, n)))

  w <- runif(n, 0.5, 2)
  w <- w / sum(w)
  Qn <- torch::torch_tensor(w)
  expect_tensor_close(Q_fluctuation(eps0, Kmat, Qn), torch::torch_tensor(w))
  expect_equal(as.numeric(Q_fluctuation(eps0, Kmat, Qn)$sum()), 1, tolerance = 1e-6)
})


test_that("dQ_fluctuation_depsilon is the Jacobian of Q_fluctuation", {
  n <- 8L; p <- 2L
  set.seed(6)
  Kmat <- matrix(rnorm(n * p), n, p)
  w <- runif(n, 0.5, 2)
  w <- w / sum(w)
  eps <- c(0.11, -0.07)

  Kt <- torch::torch_tensor(Kmat, dtype = torch::torch_double())
  Qt <- torch::torch_tensor(w, dtype = torch::torch_double())

  got <- as.matrix(
    dQ_fluctuation_depsilon(torch::torch_tensor(eps, dtype = torch::torch_double()), Kt, Qt)$detach()
  )

  expect_lt(max(abs(colSums(got))), 1e-12)

  expect_equal(
    got,
    autograd_jacobian(function(e) Q_fluctuation(e, Kt, Qt), eps, dtype = torch::torch_double()),
    tolerance = 1e-10
  )

  f_dbl <- function(e) {
    as.numeric(Q_fluctuation(torch::torch_tensor(e, dtype = torch::torch_double()), Kt, Qt))
  }

  expect_equal(got, numeric_jacobian(f_dbl, eps), tolerance = 1e-6)

  # Check the float32 path is consistent at float32 tolerance
  Kf <- torch::torch_tensor(Kmat)
  Qf <- torch::torch_tensor(w)
  expect_equal(
    as.matrix(dQ_fluctuation_depsilon(torch::torch_tensor(eps), Kf, Qf)$detach()),
    got, tolerance = 1e-5, ignore_attr = TRUE
  )
})

test_that("B() reaches the first-order condition E_Q[Ldot(psi, beta)] = 0", {
  n <- 40L; K <- 1L; p <- 2L
  set.seed(7)
  V <- rnorm(n)
  Xa <- array(cbind(1, V), dim = c(n, K, p))
  psi_a <- matrix(0.3 + 0.5 * V + rnorm(n, sd = 0.1), n, K)
  Qv <- rep(1 / n, n)
  Lm_fn <- Lm(loss_squared_error, working_model_linear)

  beta <- B(Lm_fn, torch::torch_tensor(psi_a), torch::torch_tensor(Xa), torch::torch_tensor(Qv), p)

  g <- dL(Lm_fn, torch::torch_tensor(psi_a), leaf(as.numeric(beta)),
          torch::torch_tensor(Xa), weight = torch::torch_tensor(Qv))[[1]]

  expect_lt(max(abs(as.numeric(g))), 1e-5)

  # For linear + squared error the minimiser is weighted OLS in closed form.
  Xm <- cbind(1, V)
  expect_equal(as.numeric(beta), as.numeric(solve(crossprod(Xm), crossprod(Xm, psi_a[, 1]))), tolerance = 1e-5)
})

test_that("normalizing matrix input is a symmetric Hessian", {
  # ddL returns E_Q[Lddot], a Hessian, and therefore symmetric.
  n <- 40L; K <- 1L; p <- 2L
  set.seed(7)
  V <- rnorm(n)
  Xa <- array(cbind(1, V), dim = c(n, K, p))
  psi_a <- matrix(0.3 + 0.5 * V + rnorm(n, sd = 0.1), n, K)
  b <- rnorm(p)
  Qv <- rep(1 / n, n)

  M <- ddL(Lm(loss_squared_error, working_model_linear),
           leaf(psi_a), leaf(b), torch::torch_tensor(Xa), p, torch::torch_tensor(Qv))
  Mm <- as.matrix(M$detach())
  expect_lt(max(abs(Mm - t(Mm))), 1e-5 * max(abs(Mm)))
})

test_that("clever directions are consistent with dB_dpsi", {
  data <- sim1_data(n = 100L)
  nu <- oracle_nuisance_dose_response(data)

  it <- dose_response_internals(data, nu, formula = ~A, tmle = FALSE)
  problem <- it$problem
  n <- problem$n
  p <- problem$p

  check <- function(Qv, label) {
    state <- it$state0
    state$Q <- torch::torch_tensor(Qv)
    psi <- it$spec$psi_from_state(problem, state)
    Minv <- normalizing_matrix(problem$Lm_fn, psi, state$beta, problem$design_matrix, state$Q, p)
    J <- dB_dpsi(problem$Lm_fn, psi, state$beta, problem$design_matrix, state$Q, p)
    d <- clever_directions(problem, psi, state$beta, Minv)
    expect_equal(
      as.array(J$permute(c(2, 1, 3))$detach()),
      as.array((-d * state$Q$reshape(c(n, 1, 1)))$detach()),
      tolerance = 1e-5, info = label
    )
  }

  check(rep(1 / n, n), "uniform Q")
  w <- runif(n, 0.5, 2)
  check(w / sum(w), "non-uniform Q")
})

test_that("B_wls agrees with the LBFGS solve for CATE", {
  data <- sim1_data(n = 100L)
  nu <- oracle_nuisance_cate(data)
  it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = "linear"))

  pre <- new_B_wls(it$problem)

  expect_false(is.null(pre))

  for(i in seq_len(10L)) {
    Q <- runif(it$problem$n, 0.5, 2)
    Q <- Q / sum(Q)

    psi <- matrix(rnorm(it$problem$n * it$problem$K), it$problem$n, it$problem$K)
    fast <- B_wls(pre, psi, Q)$beta
    slow <- as.numeric(B(it$problem$Lm_fn, as_float_tensor(psi)$requires_grad_(TRUE), it$problem$design_matrix, torch::torch_tensor(Q), it$problem$p))
    expect_equal(fast, slow, tolerance = 1e-4)

    g <- dL(it$problem$Lm_fn, as_float_tensor(psi), leaf(fast), it$problem$design_matrix, weight = torch::torch_tensor(Q))[[1]]
    expect_lt(max(abs(as.numeric(g))), 1e-6)
  }
})

test_that("dQ_deps_wls matches dQ_fluctuation_depsilon", {
  set.seed(1)
  n <- 40L
  p <- 2L
  Kf <- matrix(rnorm(n * p), n, p)
  Q <- runif(n, 0.5, 2)
  Q <- Q / sum(Q)
  for(e in list(rep(0, p), c(0.11, -0.08), c(-0.4, 0.5))) {
    Qt <- Q_fluctuation(torch::torch_tensor(e), torch::torch_tensor(Kf), torch::torch_tensor(Q))
    got <- dQ_deps_wls(Kf, as.numeric(Qt))
    want <- as.matrix(dQ_fluctuation_depsilon(torch::torch_tensor(e), torch::torch_tensor(Kf), torch::torch_tensor(Q))$detach())
    expect_equal(got, want, tolerance = 1e-5, ignore_attr = TRUE)
    expect_lt(max(abs(colSums(got))), 1e-6)
  }
})

test_that("bayes_jacobian_wls matches the torch Jacobian", {
  data <- sim1_data(n = 100L)
  nu <- oracle_nuisance_cate(data)
  it <- cate_internals(data, nu, formula = ~X4, tmle = TRUE)
  problem <- it$problem
  spec <- it$spec

  pre <- problem$B_wls
  pre$K_Q <- as.matrix(it$fit$final$K_Q)
  step <- spec$steps(problem)[[1]]
  step$fluctuate_Q <- TRUE

  clever <- scale_bayes_clever(problem, spec, it$fit$state, it$fit$final$clever, linear = TRUE)

  for(e in list(c(0.02, -0.03), c(-0.05, 0.03), c(0.08, 0.06))) {
    et <- torch::torch_tensor(e)
    st <- spec$apply_update(problem, step, it$fit$state, et, clever, it$fit$final$K_Q, linear = TRUE)

    psi <- spec$psi_from_state(problem, st)
    pm <- as.matrix(psi$detach())
    sol <- B_wls(pre, pm, as.numeric(st$Q))
    b <- torch::torch_tensor(sol$beta, requires_grad = TRUE)
    dps <- spec$dpsi_depsilon(problem, st, clever, et, linear = TRUE)

    slow <- bayes_jacobian(problem, psi, it$fit$state$Q, st$Q, b, dps, et, it$fit$final$K_Q)
    fast <- bayes_jacobian_wls(pre, sol, pm, st$Q, dps, dQ_deps_wls(pre$K_Q, st$Q))

    expect_equal(fast, as.matrix(slow$detach()), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(log_abs_det(fast), log_abs_det(slow), tolerance = 1e-5)
  }
})

# ----- EIF -----

test_that("eif() equals the closed-form EIF elementwise (~1)", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)
  n <- nrow(data)

  psi <- nu$mu1 - nu$mu0
  Delta <- cate_delta(data, nu)
  Xmat <- model.matrix(~ 1, data = data)

  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~1, tmle = FALSE, bayes = FALSE, nuisance_estimates = nu)

  beta <- res$plugin$est

  expect_equal(as.matrix(res$onestep$eif),
               analytic_eif_linear(Xmat, psi, beta, Delta),
               tolerance = 1e-5, ignore_attr = TRUE)
})


test_that("eif() equals the closed-form EIF elementwise (~X4)", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)
  n <- nrow(data)

  psi <- nu$mu1 - nu$mu0
  Delta <- cate_delta(data, nu)
  Xmat <- model.matrix(~ X4, data = data)

  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~1 + X4, tmle = FALSE, bayes = FALSE, nuisance_estimates = nu)

  beta <- res$plugin$est

  expect_equal(as.matrix(res$onestep$eif),
               analytic_eif_linear(Xmat, psi, beta, Delta),
               tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("cate TMLE solves the EIF estimating equation", {
  # Verify Pn[D*(mu*, Q*, eta)] \approx 0

  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)

  for(linear in c(TRUE, FALSE)) {
    res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
                formula = ~X4, tmle = tmle_control(fluctuation = if(linear) "linear" else "logistic"),
                outcome = "binomial",
                nuisance_estimates = nu)

    expect_solves_eif(res$tmle$eif)
  }
})

test_that("dose_response TMLE solves the EIF estimating equation", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_dose_response(data)

  res <- dose_response(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~A, tmle = TRUE,
              outcome = "binomial",
              bayes = FALSE, nuisance_estimates = nu)

  expect_solves_eif(res$tmle$eif)
})

# ----- ATE reduction -----
# With V = empty (formula = ~1) and squared-error loss, the CATE NP-MSM
# reduces exactly to the ATE: B(P) = E_P[mu1(X) - mu0(X)]

test_that("plug-in with ~1 equals the plug-in ATE", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)

  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~1, tmle = FALSE, bayes = FALSE, nuisance_estimates = nu)

  expect_equal(as.numeric(res$plugin$est), mean(nu$mu1 - nu$mu0), tolerance = 1e-6)
})


test_that("one-step with ~1 equals the AIPW ATE", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)

  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~1, tmle = FALSE, bayes = FALSE, nuisance_estimates = nu)

  psi <- nu$mu1 - nu$mu0
  Delta <- cate_delta(data, nu)
  expect_equal(as.numeric(res$onestep$est), mean(psi) + mean(Delta), tolerance = 1e-5)
  expect_equal(as.numeric(res$onestep$se), stats::sd(psi - mean(psi) + Delta) / sqrt(nrow(data)), tolerance = 1e-4)
})


test_that("TMLE with ~1 agrees with the AIPW ATE", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)

  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~1, tmle = TRUE, bayes = FALSE, nuisance_estimates = nu)

  aipw <- mean(nu$mu1 - nu$mu0) + mean(cate_delta(data, nu))
  expect_lt(abs(as.numeric(res$tmle$est) - aipw), 3 * as.numeric(res$tmle$se))
  expect_solves_eif(res$tmle$eif)
})

# ----- Design decisions -----
test_that("working models with p != d are supported", {
  working_model_power_law <- function(beta, X) {
    beta[1] * torch::torch_pow(X[, , 1], beta[2])
  }

  data <- sim1_data(n = 150L)
  data$f <- abs(data$X4) + 1
  nu <- oracle_nuisance_cate(data)
  res <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y",
              formula = ~-1 + f, working_model = working_model_power_law,
              p = 2L, tmle = FALSE, bayes = FALSE, nuisance_estimates = nu)

  expect_equal(res$p, 2L)
  expect_equal(res$d, 1L)
  expect_length(res$plugin$est, 2L)
  expect_equal(dim(as.matrix(res$onestep$eif)), c(nrow(data), 2L))
})

test_that("every mu component responds to epsilon, and eps = 0 is identity", {
  for(linear in c(TRUE, FALSE)) {
    for(param in c("cate", "dose_response")) {
      data <- sim1_data(n = 60L, seed = 10016)
      if(param == "cate") {
        nu <- oracle_nuisance_cate(data)
        it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = if(linear) "linear" else "logistic"))
      } else {
        nu <- oracle_nuisance_dose_response(data)
        it <- dose_response_internals(data, nu, formula = ~A, tmle = tmle_control(fluctuation = if(linear) "linear" else "logistic"))
      }

      problem <- it$problem
      spec <- it$spec
      step <- spec$steps(problem)[[1]]
      state0 <- it$state0

      # (M1a): fluctuation is identity at epsilon = 0

      z <- spec$apply_update(problem, step, state0, torch::torch_zeros(problem$p), it$clever0, it$K_Q0, linear = linear)

      for(name in names(state0$mu)) {
        expect_equal(as.array(z$mu[[name]]), as.array(state0$mu[[name]]),
                     tolerance = 1e-6, info = name)
      }

      # Every component must actually move.
      eps <- torch::torch_tensor(rep(0.25, problem$p))
      s1 <- spec$apply_update(problem, step, state0, eps, it$clever0, it$K_Q0, linear = linear)
      for(name in names(state0$mu)) {
        expect_gt(max(abs(as.array(s1$mu[[name]]) - as.array(state0$mu[[name]]))), 1e-1) # mu should move a lot
      }
      expect_gt(max(abs(as.numeric(s1$Q) - as.numeric(state0$Q))), 1e-3) # Q might not fluctuate much
    }
  }
})

# ----- generalized Bayes -----

test_that("The Jacobian dbeta/deps equals P0[lambda*]", {
  data <- sim1_data(n = 500L)
  nu <- oracle_nuisance_cate(data)
  it <- cate_internals(data, nu, formula = ~X4, tmle = tmle_control(fluctuation = "logistic"))
  n <- it$problem$n
  p <- it$problem$p
  fn <- it$fit$final

  e0 <- torch::torch_tensor(rep(0, p))
  st <- it$spec$apply_update(it$problem, list(id = 1L, fluctuate_Q = TRUE),
                             it$fit$state, e0, fn$clever, fn$K_Q, linear = FALSE)

  psi <- it$spec$psi_from_state(it$problem, st)
  beta <- B(it$problem$Lm_fn, psi, it$problem$design_matrix, st$Q, p)
  dpsi <- it$spec$dpsi_depsilon(it$problem, st, fn$clever, e0, linear = FALSE)
  J <- as.matrix(bayes_jacobian(it$problem, psi, it$fit$state$Q, st$Q, beta, dpsi, e0, fn$K_Q)$detach())

  # Analytical jacobian for CATE + linear working model + squared error loss
  X   <- cbind(1, data$X4)
  Sig <- crossprod(X) / n
  Si <- solve(Sig)
  m1  <- as.numeric(st$mu$a1)
  m0 <- as.numeric(st$mu$a0)
  pp <- nu$pi
  w   <- m1 * (1 - m1) / pp + m0 * (1 - m0) / (1 - pp)
  r   <- (m1 - m0) - as.vector(X %*% it$tmle$est)

  Jpsi_ref <- -Si %*% (crossprod(X * sqrt(w))   / n) %*% Si
  JQ_ref   <- -Si %*% (crossprod(X * abs(r))    / n) %*% Si

  # Same split from the code, using the include_Q switch
  e0   <- torch::torch_tensor(rep(0, it$problem$p))
  psi  <- it$spec$psi_from_state(it$problem, st)
  beta <- B(it$problem$Lm_fn, psi, it$problem$design_matrix, st$Q, it$problem$p)
  dpsi <- it$spec$dpsi_depsilon(it$problem, st, fn$clever, e0, linear = FALSE)

  Jpsi <- as.matrix(bayes_jacobian(it$problem, psi, st$Q, st$Q, beta, dpsi, e0,
                                   fn$K_Q, include_Q = FALSE)$detach())
  Jall <- as.matrix(bayes_jacobian(it$problem, psi, st$Q, st$Q, beta, dpsi, e0,
                                   fn$K_Q, include_Q = TRUE)$detach())

  expect_equal(J, Jpsi_ref + JQ_ref, tolerance = 5e-3, ignore_attr = TRUE)
})
