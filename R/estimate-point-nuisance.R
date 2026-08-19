#' Estimate nuisance parameters for point-treatment NP-MSMs
#'
#' Cross-fits the propensity score(s) and the outcome regression for a
#' single-time-point target functional. The CATE (`K = 2`) and the categorical
#' dose-response (`K > 2`) functionals are both special cases.
#'
#' @noRd
estimate_point_nuisance <- function(
  data,
  X,
  A,
  Y,
  control,
  outcome_type,
  levels = NULL,
  propensity = c("one_vs_rest", "binary"),
  estimate_conditional_variance = FALSE
) {
  propensity <- match.arg(propensity)
  n <- nrow(data)

  setup <- nuisance_setup(n, control, outcome_type)

  if(is.null(levels)) levels <- sort(unique(data[[A]]))
  K <- length(levels)

  if(propensity == "binary" && K != 2L) {
    stop('propensity = "binary" requires exactly two treatment levels.')
  }

  if(!all(data[[A]] %in% levels)) {
    stop("Some observed values of `", A, "` are not present in `levels`.")
  }

  pi_a_hat <- mu_a_hat <- condvar_a_hat <- matrix(0, n, K)

  for (fold in seq_along(setup$cv)) {
    train <- setup$cv[[fold]]$training_set
    valid <- setup$cv[[fold]]$validation_set

    # ----- Propensity score(s) -----
    if(propensity == "binary") {
      res <- sl_fit_predict(
        Y = as.numeric(data[[A]][train] == levels[2]),
        X = data[train, X, drop = FALSE],
        newdata = list(valid = data[valid, c(X), drop = FALSE]),
        SL.library = control$learners_trt,
        family = "binomial",
        cv_control = setup$cv_control,
        bounds = c(0, 1),
        epsilon = control$epsilon
      )

      pi_a_hat[valid, 2] <- res$pred$valid
      pi_a_hat[valid, 1] <- 1 - res$pred$valid
    }
    else {
      for (k in seq_len(K)) {
        res <- sl_fit_predict(
          Y = as.numeric(data[[A]][train] == levels[2]),
          X = data[train, X, drop = FALSE],
          newdata = list(valid = data[valid, c(X), drop = FALSE]),
          SL.library = control$learners_trt,
          family = "binomial",
          cv_control = setup$cv_control,
          bounds = c(0, 1),
          epsilon = control$epsilon
        )

        pi_a_hat[valid, k] <- res$pred$valid
      }
    }

    # ----- Outcome regresion -----
    newdata <- lapply(seq_len(K), function(k) {
      newdata <- data[valid, c(X, A), drop = FALSE]
      newdata[[A]] <- levels[k]
      newdata
    })
    newdata[["train"]] <- data[train, c(X, A), drop = FALSE]

    res <- sl_fit_predict(
      Y = data[[Y]][train],
      X = data[train, c(X, A), drop = FALSE],
      newdata = newdata,
      SL.library = control$learners_outcome,
      family = setup$outcome_family,
      bounds = if(outcome_type == "binomial") c(0, 1) else NULL,
      cv_control = setup$cv_control
    )

    for (k in seq_len(K)) {
      mu_a_hat[valid, k] <- res$pred[[k]]
    }

    if(identical(outcome_type, "binomial")) {
      # For binary outcomes, closed form for conditional variance is available
      for(k in seq_len(K)) {
        condvar_a_hat[valid, k] <- mu_a_hat[valid, k] * (1 - mu_a_hat[valid, k])
      }
    }
    else if (isTRUE(estimate_conditional_variance)) {
      # Otherwise, only estimate conditional variance if requested
      yhat <- res$pred[["train"]]
      y2 <- (data[[Y]][train] - yhat)^2

      yvar <- sl_fit_predict(
        Y = y2,
        X = data[train, c(X, A), drop = FALSE],
        newdata = newdata,
        SL.library = control$learners_outcome,
        family = setup$outcome_family,
        cv_control = setup$cv_control
      )

      for(k in seq_len(K)) {
        condvar_a_hat[valid, k] <- ifelse(yvar$pred[[k]] < control$epsilon, control$epsilon, yvar$pred[[k]])
      }
    }
  }

  # ------ Observed-arm nuisances -----
  mu_hat <- pi_obs_hat <- condvar_hat <- numeric(n)
  for (k in seq_len(K)) {
    ind <- data[[A]] == levels[k]
    mu_hat[ind] <- mu_a_hat[ind, k]
    condvar_hat[ind] <- condvar_a_hat[ind, k]
    pi_obs_hat[ind] <- pi_a_hat[ind, k]
  }

  list(
    levels = levels,
    K = K,
    cv = setup$cv,
    pi_a = pi_a_hat,       # n x K: P(A = levels[k] | X)
    pi_obs = pi_obs_hat,   # x:     P(A = A_i | X)
    mu_a = mu_a_hat,       # n x K: E[Y | A = levels[k], X]
    mu = mu_hat,           # n:     E[Y | A = A_i, X]
    condvar_a = condvar_a_hat,
    condvar = condvar_hat
  )
}
