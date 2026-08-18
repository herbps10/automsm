#' Group identical columns of a matrix
#'
#' At node t, the regression depends on the target column alone,
#' because A_t only matters for *prediction*. Two regimes
#' needing the same fit are those with identical targets.
#'
#'   t == tau: every column is Y      -> 1 fit, K predictions
#'   t < tau:  untargeted chain       -> 1 fit per distinct suffix
#'   t < tau:  targeted (TMLE) chain  -> K fits
#'
#' @noRd
column_keys <- function(mat) {
  K <- ncol(mat)
  key <- integer(K)
  nxt <- 0L
  for(j in seq_len(K)) {
    for(g in seq_len(j - 1L)) {
      if(identical(mat[, j], mat[, g])) {
        key[j] <- key[g]
        break
      }
    }
    if(key[j] == 0L) {
      nxt <- nxt + 1L
      key[j] <- nxt
    }
  }
  key
}

#' Fit one sequential regression within a single fold
#'
#' Regress \code{target_train} on (H_t, A_t) using the training rows of \code{fold},
#' then evaluate at A_t = a_t^(j) for every regime j, on both training and validation folds.
#' The next node down regresses on the training predictions, and the pooled fluctuations
#' use the validation predictions.
#'
#' @param target_train |T_v| x K; column j is regime j's current pseudo-outcome
#'   on the training rows. At t = tau this is Y in every column.
#' @noRd
fit_ice_node <- function(t, target_train, data, Ls, As, regimes, fold,
                         learners, family, bounds, epsilon, cv_control) {
  train <- fold$training_set
  valid <- fold$validation_set

  K <- nrow(regimes)
  stopifnot(nrow(target_train) == length(train), ncol(target_train) == K)

  Ht_At <- c(history_cols(t, Ls, As), As[t])
  Xtrain <- data[train, Ht_At, drop = FALSE]

  key <- column_keys(target_train)
  out_train <- matrix(NA_real_, length(train), K)
  out_valid <- matrix(NA_real_, length(valid), K)
  fits <- vector("list", max(key))

  for(g in seq_len(max(key))) {
    j0 <- which(key == g)[1L]
    fits[[g]] <- sl_fit(Y = target_train[, j0], X = Xtrain,
                        SL.library = learners, family = family,
                        cv_control = cv_control)

    for(j in which(key == g)) {
      nd_train <- data[train, Ht_At, drop = FALSE]
      nd_valid <- data[valid, Ht_At, drop = FALSE]
      nd_train[[As[t]]] <- regimes[j, As[t]]
      nd_valid[[As[t]]] <- regimes[j, As[t]]
      out_train[, j] <- sl_predict(fits[[g]], nd_train, bounds, epsilon)
      out_valid[, j] <- sl_predict(fits[[g]], nd_valid, bounds, epsilon)
    }
  }

  list(train = out_train, valid = out_valid, fits = fits, key = key)
}

#' Untargeted sequential regression chain
#'
#' outer loop over folds, inner fold goes backwward over t,
#' with each node's target being the previous (t + 1) node's training predictions.
#' @noRd
estimate_ice_chain <- function(data, Ls, As, Y, regimes, folds, learners, outcome_family,
                               cv_control, bounds, epsilon) {
  n <- nrow(data)
  tau <- length(As)
  K <- nrow(regimes)
  mu <- array(NA_real_, dim = c(K, n, tau + 1L))

  for(j in seq_len(K)) {
    mu[j, , tau + 1L] <- data[[Y]]
  }

  for(v in seq_along(folds)) {
    train <- folds[[v]]$training_set
    valid <- folds[[v]]$validation_set
    target <- matrix(data[[Y]][train], nrow = length(train), ncol = K)
    for(t in tau:1) {
      family <- if(t == tau) outcome_family else stats::gaussian()
      nd <- fit_ice_node(t, target, data, Ls, As, regimes, folds[[v]],
                         learners, family, bounds, epsilon, cv_control)

      mu[, valid, t] <- t(nd$valid)
      target <- nd$train
    }
  }
  stopifnot(!anyNA(mu))
  mu
}
