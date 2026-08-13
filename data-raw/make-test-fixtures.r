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
  nuisance <- estimate_cate_nuisance(
    data = dat,
    X = c("X1", "X2"),
    A = "A",
    Y = "Y",
    learners_trt = "SL.glm.interaction",
    learners_outcome = "SL.glm.interaction",
    outer_folds = 2,
    inner_folds = 2,
    outcome_type = outcome_type,
    estimate_conditional_variance = (outcome_type == "continuous"),
    epsilon = 1e-5
  )

  fit_onestep <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = FALSE,
    nuisance = nuisance
  )

  fit_tmle <- cate(
    dat,
    c("X1", "X2"),
    "A",
    "Y",
    formula = ~1 + X2,
    outcome_type = outcome_type,
    tmle = TRUE,
    nuisance = nuisance
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
    bayes = TRUE,
    bayes_chains = 2,
    bayes_draws = 100,
    nuisance = nuisance
  )

  saveRDS(
    list(
      dat = dat,
      nuisance = nuisance,
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
  nuisance <- estimate_longitudinal_dose_response_nuisance(
    data = dat,
    Ls = Ls,
    As = As,
    Y = "Y",
    regimes = regimes,
    learners_trt = "SL.glm.interaction",
    learners_outcome = "SL.glm.interaction",
    outer_folds = 2,
    inner_folds = 2,
    outcome_type = outcome_type,
    epsilon = 1e-5,
    cv = NULL
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
    nuisance = nuisance
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
    nuisance = nuisance
  )

  saveRDS(
    list(
      dat = dat,
      nuisance = nuisance,
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

