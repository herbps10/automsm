# Shared functions for single-time point MSMs (cate and dose_response)



#' Estimate nuisance parameters for point-treatment NP-MSMs
#'
#' Cross-fits the propensity score(s) and the outcome regression for a
#' single-time-point target functional. The CATE (`K = 2`) and the categorical
#' dose-response (`K > 2`) functionals are both special cases.
#'
#' @param data A `data.frame`.
#' @param X Character vector of covariate column names.
#' @param A String naming the treatment column.
#' @param Y String naming the outcome column.
#' @param learners_trt SuperLearner library for the propensity score(s).
#' @param learners_outcome SuperLearner library for the outcome regression.
#' @param outer_folds Number of cross-fitting folds.
#' @param inner_folds Number of inner SuperLearner CV folds.
#' @param outcome_type Either \code{"continuous"} or \code{"binomial"}.
#' @param levels Treatment levels. Defaults to \code{sort(unique(data[[A]]))}.
#' @param propensity Either \code{"binary"} (fit a single model for \code{levels[2]}
#'   and derive \code{levels[1]} as its complement) or \code{"one_vs_rest"} (fit a
#'   separate model per level). Note these are not numerically equivalent, even
#'   when \code{K = 2}.
#' @param estimate_conditional_variance Whether to fit the conditional variance.
#' @param epsilon Bounding constant for propensities/binomial means.
#' @param cv Optional pre-computed \code{origami} fold object, to share
#'   folds across estimators.
#' @noRd
estimate_point_nuisance <- function(
  data,
  X,
  A,
  Y,
  learners_trt,
  learners_outcome,
  outer_folds,
  inner_folds,
  outcome_type,
  levels = NULL,
  propensity = c("one_vs_rest", "binary"),
  estimate_conditional_variance = FALSE,
  epsilon = 1e-5,
  cv = NULL
) {
  propensity <- match.arg(propensity)
  n <- nrow(data)

  setup <- nuisance_setup(n, outer_folds, inner_folds, outcome_type, cv)

  if(is.null(levels)) levels <- sort(unique(data[[A]]))
  K <- length(levels)

  if(propensity == "binary" && K != 2L) {
    stop('propensity = "binary" requires exactly two treatment levels.')
  }

  if(!all(data[[A]] %in% levels)) {
    stop("Some observed values of `", A, "` are not present in `levels`.")
  }

  pi_a_hat <- mu_a_hat <- matrix(0, n, K)
  condvar_hat <- numeric(n)

  for (fold in seq_along(setup$cv)) {
    train <- setup$cv[[fold]]$training_set
    valid <- setup$cv[[fold]]$validation_set

    # ----- Propensity score(s) -----
    if(propensity == "binary") {
      res <- sl_fit_predict(
        Y = as.numeric(data[[A]][train] == levels[2]),
        X = data[train, X, drop = FALSE],
        newdata = list(valid = data[valid, c(X), drop = FALSE]),
        SL.library = learners_trt,
        family = "binomial",
        cv_control = setup$cv_control,
        bounds = c(0, 1),
        epsilon = epsilon
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
          SL.library = learners_trt,
          family = "binomial",
          cv_control = setup$cv_control,
          bounds = c(0, 1),
          epsilon = epsilon
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
      SL.library = learners_outcome,
      family = setup$outcome_family,
      bounds = if(outcome_type == "binomial") c(0, 1) else NULL,
      cv_control = setup$cv_control
    )

    for (k in seq_len(K)) {
      mu_a_hat[valid, k] <- res$pred[[k]]
    }

    if (isTRUE(estimate_conditional_variance)) {
      yhat <- res$pred[["train"]]
      y2 <- (data[[Y]][train] - yhat)^2

      yvar <- sl_fit_predict(
        Y = y2,
        X = data[train, c(X, A), drop = FALSE],
        newdata = list(valid = data[valid, c(X, A), drop = FALSE]),
        SL.library = learners_outcome,
        family = setup$outcome_family,
        cv_control = setup$cv_control
      )

      condvar_hat[valid] <- ifelse(yvar$pred$valid < epsilon, epsilon, yvar$pred$valid)
    }
  }

  # ------ Observed-arm nuisances -----
  mu_hat <- pi_obs_hat <- numeric(n)
  for (k in seq_len(K)) {
    ind <- data[[A]] == levels[k]
    mu_hat[ind] <- mu_a_hat[ind, k]
    pi_obs_hat[ind] <- pi_a_hat[ind, k]
  }

  list(
    levels = levels,
    K = K,
    cv = cv,
    pi_a = pi_a_hat,       # n x K: P(A = levels[k] | X)
    pi_obs = pi_obs_hat,   # x:     P(A = A_i | X)
    mu_a = mu_a_hat,       # n x K: E[Y | A = levels[k], X]
    mu = mu_hat,           # n:     E[Y | A = A_i, X]
    condvar = condvar_hat
  )
}
