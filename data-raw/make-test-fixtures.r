make_cate_fixture <- function(path, binary, nonlinear) {
  dat <- simulate_cate(
    N = 500,
    sigma = 0.1,
    seed = 1,
    nonlinear = nonlinear,
    binary = binary
  )

  outcome_type <- if(isTRUE(binary)) "binomial" else "continuous"

  set.seed(1)
  nuisance_estimates <- estimate_cate_nuisance(
    data = dat,
    X = c("X1", "X2"),
    A = "A",
    Y = "Y",
    control = nuisance_control(
      learners_trt = "SL.glm.interaction",
      learners_outcome = "SL.glm.interaction",
      outer_folds = 2,
      inner_folds = 2,
      epsilon = 1e-5
    ),
    outcome_type = outcome_type,
    estimate_conditional_variance = (outcome_type == "continuous")
  )

  fit_onestep <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = FALSE,
    nuisance_estimates = nuisance_estimates
  )

  fit_tmle <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = TRUE,
    nuisance_estimates = nuisance_estimates
  )

  set.seed(1)
  fit_bayes <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100,
      warmup = 100
    ),
    nuisance_estimates = nuisance_estimates
  )

  saveRDS(
    list(
      dat = dat,
      nuisance_estimates = nuisance_estimates,
      fit_onestep = fit_onestep,
      fit_tmle = fit_tmle,
      fit_bayes = fit_bayes
    ),
    path
  )
}

make_dose_response_fixture <- function(path, binary, nonlinear) {
  dat <- simulate_dose_response(
    N = 500,
    sigma = 0.1,
    seed = 1,
    treatments = 6,
    nonlinear = nonlinear,
    binary = binary
  )

  outcome_type <- if(isTRUE(binary)) "binomial" else "continuous"

  set.seed(1)
  nuisance_estimates <- estimate_dose_response_nuisance(
    data = dat,
    X = c("X1", "X2"),
    A = "A",
    Y = "Y",
    control = nuisance_control(
      learners_trt = "SL.glm.interaction",
      learners_outcome = "SL.glm.interaction",
      outer_folds = 2,
      inner_folds = 2,
      epsilon = 1e-5
    ),
    outcome_type = outcome_type,
    estimate_conditional_variance = (outcome_type == "continuous")
  )

  fit_onestep <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = FALSE,
    nuisance_estimates = nuisance_estimates
  )

  fit_tmle <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = TRUE,
    nuisance_estimates = nuisance_estimates
  )

  set.seed(1)
  fit_bayes <- dose_response(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = TRUE,
    bayes = bayes_control(
      chains = 2,
      draws = 100,
      warmup = 100
    ),
    nuisance_estimates = nuisance_estimates
  )

  saveRDS(
    list(
      dat = dat,
      nuisance_estimates = nuisance_estimates,
      fit_onestep = fit_onestep,
      fit_tmle = fit_tmle,
      fit_bayes = fit_bayes
    ),
    path
  )
}

make_longitudinal_dose_response_fixture <- function(path, binary) {
  tau <- 2
  dat <- simulate_longitudinal_dose_response(
    N = 500,
    sigma = 0.1,
    tau = tau,
    seed = 1,
    binary = binary
  )

  outcome_type <- if(isTRUE(binary)) "binomial" else "continuous"
  Ls <- lapply(seq_len(tau), function(t) paste0("L", t))
  As <- paste0("A", seq_len(tau))
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  set.seed(1)
  nuisance_estimates <- estimate_longitudinal_dose_response_nuisance(
    data = dat,
    Ls = Ls,
    As = As,
    Y = "Y",
    regimes = regimes,
    control = nuisance_control(
      learners_trt = "SL.glm.interaction",
      learners_outcome = "SL.glm.interaction",
      outer_folds = 2,
      inner_folds = 2,
      epsilon = 1e-5,
    ),
    outcome_type = outcome_type
  )

  loss <- if(outcome_type == "continuous") loss_squared_error else loss_cross_entropy_logit

  fit_onestep <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss,
    outcome_type = outcome_type,
    tmle = FALSE,
    nuisance_estimates = nuisance_estimates
  )

  set.seed(1)
  fit_tmle <- longitudinal_dose_response(
    dat,
    Ls,
    As,
    "Y",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss,
    outcome_type = outcome_type,
    tmle = TRUE,
    nuisance_estimates = nuisance_estimates
  )

  saveRDS(
    list(
      dat = dat,
      nuisance_estimates = nuisance_estimates,
      fit_onestep = fit_onestep,
      fit_tmle = fit_tmle,
      tau = tau,
      Ls = Ls,
      As = As,
      regimes = regimes
    ),
    path
  )
}

# ----- Save fixtures -----
make_cate_fixture(test_path("fixtures", "cate_continuous_linear.rds"), binary = FALSE, nonlinear = FALSE)
make_cate_fixture(test_path("fixtures", "cate_binary_linear.rds"), binary = TRUE, nonlinear = FALSE)

make_dose_response_fixture(test_path("fixtures", "dose_response_continuous_linear.rds"), binary = FALSE, nonlinear = FALSE)
make_dose_response_fixture(test_path("fixtures", "dose_response_binary_linear.rds"), binary = TRUE, nonlinear = FALSE)

make_longitudinal_dose_response_fixture(test_path("fixtures", "longitudinal_dose_response_continuous.rds"), binary = FALSE)
make_longitudinal_dose_response_fixture(test_path("fixtures", "longitudinal_dose_response_binary.rds"), binary = TRUE)

