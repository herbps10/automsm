#' Estimate a Non-parametric Marginal Structural Model for Longitudinal Dose Response
#'
#' Fits a non-parametric marginal structural model (NP-MSM) summarizing the mean
#' outcome under static longitudinal treatment trajectories \eqn{\bar{a} \in \{0, 1\}^\tau}.
#' The target functional is the set of iterated conditional expectations
#' \eqn{\psi_P^{(\bar{a})}(L_1)} under each static longitudinal policy, projected
#' onto a user-supplied working model via a user-supplied loss function. Plug-in,
#' one-step, and Targeted Minimum Loss-based (TMLE) estimators are computed.
#' The one-step and TMLE estimators are semi-parametric efficient and yield
#' asymptotically valid confidence intervals.
#'
#' @details
#' The observed data are assumed to have the longitudinal structure
#' \eqn{O = (L_1, A_1, L_2, A_2, \dots, L_\tau, A_\tau, Y)}, where \eqn{L_t} are
#' time-varying covariates, \eqn{A_t \in \{0, 1\}} are time-varying treatment indicators,
#' and \eqn{Y} is an end-of-study outcome. All \eqn{2^\tau} static treatment trajectories
#' are enumerated, and the mean outcome under each is identified via iterated conditional
#' expectations (ICE). THe ICE-based conditional means are then projected onto the working
#' model defined by \code{formula}, \code{working_model}, and \code{loss}.
#'
#' Nuisance parameters (the time-varying propensity scores and the sequential
#' outcome regressions) are estimated via cross-fitting, using the \pkg{SuperLearner}
#' libraries supplied in \code{learners_trt} and \code{learners_outcome}. The same
#' cross-fitting folds are shared between the propensity score and ICE regressions.
#' Sequential positivity is required for identification; the \code{epsilon} argument
#' bounds estimated propensities away from 0 and 1 to stabilize the inverse cumulative-propensity
#' weights that appear in the efficient influence function.
#'
#'
#' @param data A \code{data.frame} containing the columns referenced by \code{Ls},
#'   \code{As}, and \code{Y}.
#' @param Ls A list of length \eqn{tau}, whose \eqn{t}-th element contains the names of
#'   the time-varying covariates measured at timepoint \eqn{t}. Must have the same length
#'   as \code{As}.
#' @param As A character vector of length \eqn{\tau} giving the treatment column names,
#'   ordered by timepoint. All entries must be columns of \code{data}.
#' @param Y A string naming the outcome column in \code{data}.
#' @param formula A model formula specifying the marginal structural working-model design matrix.
#' @param outcome_type Outcome type, either \code{"continuous"} or \code{"binomial"}.
#' @param loss The marginal structural model loss function measuring the fidelity of the
#'   working model to the target functional (e.g. \code{loss_squared_error}).
#' @param working_model The marginal structural model working model (e.g. \code{working_model_linear}).
#' @param h A weight function specifying the weight assigned to each treatment
#'   trajectory in the summed-over-regimes loss. It must accept the `regimes`
#'   data frame (with one row per trajectory in \eqn{\{0, 1\}^\tau} and one column
#'   per treatment time point) and return a numeric vector of length `nrow(regimes)`
#'   giving the (non-negative) weight for each regime. The default,
#'   `function(regimes) rep(1, nrow(regimes))`, gives equal weight to all
#'   trajectories.
#'
#'   Note that, unless the working model is assumed correctly specified,
#'   the choice of weight function changes the target parameter being
#'   estimated (see Petersen et al., 2014). Also note that, when `h` is
#'   estimated from the data, the influence-curve based inference is valid
#'   for the estimand defined by the estimated weights, and does not account
#'   for uncertainty in the estimation of `h` itself.
#' @param learners_trt A character vector of \pkg{SuperLearner} libraries for estimating
#'   the time-varying propensity scores.
#' @param learners_outcome A character vector of \pkg{SuperLearner} libraries for
#'   estimating the sequential outcome regressions.
#' @param outer_folds Number of folds in the outer cross-fitting loop.
#' @param inner_folds Number of folds for the inner \pkg{SuperLearner} cross-validation
#'   within each outer cross-fitting fold.
#' @param tmle Logical; whether to run the TMLE estimator.
#' @param tmle_maxiter Maximum number of TMLE iterations.
#' @param tmle_linear Logical; whether to use a linear TMLE fluctuation model
#'   (\code{TRUE}) or a logistic fluctuation model (\code{FALSE}).
#' @param epsilon Adjustment bounding estimated propensities/means away from 0 and 1.
#' @param nuisance Optional list of pre-computed nuisance parameters. If \code{NULL}
#'   (the default), nuisance parameters are estimated internally via cross-fitting.
#'
#' @return An object of class \code{"automsm"}: a list with components
#'   \describe{
#'     \item{estimand}{Character string, \code{"longitudinal_dose_response"}.}
#'     \item{p}{Number of working-model coefficients.}
#'     \item{n}{Sample size.}
#'     \item{tau}{Number of treatment timepoints \eqn{\tau}.}
#'     \item{formula}{The working-model formula used.}
#'     \item{working_model, loss}{The working model and loss function used.}
#'     \item{terms}{Character vector of working-model design-matrix term names.}
#'     \item{learners_trt, learners_outcome}{The \pkg{SuperLearner} libraries used.}
#'     \item{nuisance}{The (estimated or supplied) nuisance parameters.}
#'     \item{regimes}{A \code{data.frame} of the \eqn{2^\tau} enumerated treatment trajectories.}
#'     \item{plugin}{A list with the plug-in piont estimate (\code{est}).}
#'     \item{onestep}{A list with the one-step point estimate (\code{est}), standard
#'       errors (\code{se}), confidence-interval bounds (\code{lower}, \code{upper}),
#'       the estimated efficient influence function (\code{eif}), the projected
#'       conditional means (\code{psi}), and joint draws (\code{joint_draws}).}
#'   }
#'
#' @references Petersen M, Schwab J, Gruber S, Blaser N, Schomaker M, van der Laan M.
#' (2014) Targeted Maximum Likelihood Estimation for Dynamic and Static Longitudinal
#' Marginal Structural Working Models. \emph{Journal of Causal Inference}, 2(2):147-185.
#' \doi{10.1515/jci-2013-0007}
#'
#' @seealso \code{\link{dose_response}} for the analogous estimator with a
#'   high-dimensional (categorical) point treatment.
#'
#' @examples
#' \dontrun{
#' # Simulate longitudinal data with tau = 2 treatment timepoints
#' data <- simulate_longitudinal_dose_response(N = 500, tau = 2)
#'
#' fit <- longitudinal_dose_response(
#'   data = data,
#'   Ls = list(c("L1"), c("L2")),
#'   As = c("A1", "A2"),
#'   Y = "Y",
#'   formula = ~1 + I(A1 + A2),
#'   outcome_type = "binomial"
#' )
#' }
#'
#' @importFrom torch torch_tensor torch_zeros torch_reshape nn_bce_with_logits_loss nn_mse_loss
#' @importFrom stats dnorm qnorm
#' @importFrom adaptMCMC MCMC
#'
#' @export
longitudinal_dose_response <- function(
  data,
  Ls,
  As,
  Y,
  formula,
  outcome_type = "binomial",
  loss = loss_squared_error,
  working_model = working_model_linear,
  h = function(regimes) rep(1, nrow(regimes)),
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm",
  outer_folds = 5,
  inner_folds = 5,
  tmle = TRUE,
  tmle_maxiter = 25,
  tmle_linear = TRUE,
  epsilon = 1e-5,
  nuisance = NULL
) {
  # ----- Argument checks -----
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)

  # Check As
  checkmate::assert_character(As, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(As, choices = names(data))

  # Check Ls
  checkmate::assert_list(Ls, len = length(As))
  for(Lt in Ls) checkmate::assert_subset(Lt, choices = names(data))

  # Check Y
  checkmate::assert_string(Y)
  checkmate::assert_choice(Y, choices = names(data))

  checkmate::assert_formula(formula)

  checkmate::assert_choice(outcome_type, choices = c("continuous", "binomial"))

  checkmate::assert_function(loss)
  checkmate::assert_function(working_model)

  checkmate::assert_character(learners_trt, min.len = 1, any.missing = FALSE)
  checkmate::assert_character(learners_outcome, min.len = 1, any.missing = FALSE)

  checkmate::assert_count(outer_folds, positive = TRUE)
  checkmate::assert_count(inner_folds, positive = TRUE)

  checkmate::assert_flag(tmle)
  checkmate::assert_count(tmle_maxiter, positive = TRUE)
  checkmate::assert_flag(tmle_linear)

  checkmate::assert_number(epsilon, lower = 0, upper = 0.5, finite = TRUE)

  checkmate::assert(
    checkmate::check_null(nuisance),
    checkmate::check_list(nuisance),
    combine = "or"
  )


  n <- nrow(data)
  tau <- length(As)
  stopifnot(length(Ls) == length(As))

  # Enumerate all binary treatment trajectories in {0,1}^\tau
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As
  k <- nrow(regimes) # k = 2^\tau

  # Evaluate the user-specified weight function over all regimes
  h_weights <- h(regimes)
  stopifnot(
    length(h_weights) == k,
    all(is.finite(h_weights)),
    all(h_weights) >= 0
  )

  h_weights_t <- torch::torch_tensor(h_weights, dtype = torch::torch_float())

  # Map each observation to the regime index matching its observed trajectory
  regime_index <- match(
    interaction(data[, As, drop = FALSE]),
    interaction(regimes)
  )

  # Set up cross-train so the same training and valid folds
  # are used in propensity scores and ICE regressions
  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  if (outer_folds == 1) {
    cv[[1]]$training_set <- cv[[1]]$validation_set
  }

  if (is.null(nuisance)) {
    nuisance <- estimate_longitudinal_dose_response_nuisance(
      data,
      Ls,
      As,
      Y,
      regimes,
      learners_trt,
      learners_outcome,
      cv,
      outer_folds,
      inner_folds,
      outcome_type,
      epsilon = epsilon
    )
  }


  # Form the 3D MSM design matrix (K x n x p)
  # One n x p design block per regime, substituting the regime's trajectory
  mat0 <- model.matrix(formula, data = data)
  terms <- colnames(mat0)
  p <- ncol(mat0)

  Yt <- torch::torch_tensor(data[[Y]], dtype = torch::torch_float())

  design_matrix <- torch::torch_zeros(c(k, n, p))
  for(j in 1:k) {
    dataj <- data
    dataj[, As] <- regimes[j, ]
    matj <- stats::model.matrix(formula, data = dataj)
    design_matrix[j, , ] <- matj
  }

  # Combined (summed over regimes) loss
  combined_loss <- function(t, beta, X) {
    n <- rev(dim(X))[2]
    K <- dim(X)[1]
    sum <- torch::torch_tensor(rep(0, n), requires_grad = TRUE)
    for(j in 1:K) {
      # regime j's contribution is weighted by h(regime_j)
      sum <- sum$add(loss(working_model(beta, X[j, , ]), t[, j]))$mul(h_weights_t[j])
    }
    sum
  }

  ##### Plugin estimator
  Q <- torch::torch_tensor(rep(1 / n, n))

  # psi is n x K: the ICE conditional means psi_P^{(\bar{a})}(L_1) per regime.
  psi <- torch::torch_tensor(t(nuisance$mu[, , 1]), requires_grad = TRUE)

  plugin <- B(combined_loss, psi, design_matrix, Q)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  # Delta: longitudinal ICE-based EIF term (n x K)
  # Telescoping sum over t of (mu_{t+1} - mu_t) / g_{0:t}
  pi_cumprod <- cumulative_propensity_scores(regimes, nuisance$pi) # k x n x tau

  A_obs <- as.matrix(data[, As, drop = FALSE]) # k x n x tau
  stopifnot(identical(colnames(regimes), As))

  W <- on_protocol_weights(regimes, A_obs)

  resid <- (nuisance$mu[, , 2:(tau + 1), drop = FALSE] - nuisance$mu[, , 1:tau, drop = FALSE]) / pi_cumprod
  resid <- resid * W # apply on-protocol gating

  Delta <- apply(resid, c(1, 2), sum)
  Delta <- t(Delta)


  stopifnot(all(is.finite(Delta)))

  ##### One-step estimator
  onestep_est <- onestep(
    combined_loss,
    psi,
    beta,
    design_matrix,
    Q,
    torch::torch_tensor(Delta)
  )

  draws <- 1e3
  onestep_joint <- mvtnorm::rmvnorm(
    draws,
    mean = as.numeric(onestep_est$est),
    sigma = var(as.matrix(onestep_est$eif)) / n
  )

  if (tmle == TRUE) {

    # Inverse-weight clever-covariate multipliers per regime, node:
    #   HA_node[j, i, t] = I(on-protocol through t) / g_{0:t}
    HA_node <- W / pi_cumprod
    stopifnot(all(is.finite(HA_node)))

    # Clever covariate at a single ICE node k, for the current beta/psi
    calculate_clever_node <- function(Lm, HA_node_k, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]

      cleverA <- torch::torch_tensor(array(0, dim = c(k, n, p)))
      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      for (i in 1:n) {
        d <- Minv$matmul(grad_dL(
          Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE]
        ))
        for (j in 1:k) {
          x <- torch::torch_zeros(k)
          x[j] <- HA_node_k[j, i] # regime j's inverse-weight at this node
          cleverA[j, i, ] <- torch::torch_reshape(d$matmul(x), p)
        }
      }
      cleverA   # k x n x p
    }

    # Clever covariate for marginal distribution of covariates
    calculate_K <- function(Lm, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]
      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      Kmat <- torch::torch_tensor(matrix(0, n, p))
      for (i in 1:n) {
        Kmat[i, ] <- Minv$matmul(dL(
          Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE]
        )[[1]])$reshape(p)
      }
      Kmat
    }

    Q_fluctuation <- function(epsilon, K, Q) {
      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      exp(K$matmul(epsilon) * Q) / Qn
    }

    if (tmle_linear == TRUE) {
      tmle_loss <- nn_mse_loss(reduction = "sum")
    } else {
      tmle_loss <- nn_bce_with_logits_loss(reduction = "sum")
    }

    # Fluctuation object for a SINGLE node: regress the forward ICE fit
    # (mu_next, per regime) on the node's clever covariate, offset by
    # current fit.
    node_fluctuation_model <- function(epsilon, mu_node, mu_next, cleverA_node, clever_K, Q) {
      target <- 0
      for(j in 1:k) {
        if (tmle_linear == TRUE) {
          pred_j <- mu_node[, j] + cleverA_node[j, , ]$matmul(epsilon)
        } else {
          pred_j <- mu_node[, j]$logit() + cleverA_node[j, , ]$matmul(epsilon)
        }

        target <- target + tmle_loss(pred_j, mu_next[, j])
      }

      Qf <- Q_fluctuation(epsilon, clever_K, Q)
      target - log(Q)$sum()
    }

    # Initialize TMLE
    Q_star <- Q
    mu_star_nodes <- vector("list", tau + 1)
    for(t in 1:(tau + 1)) {
      mu_star_nodes[[t]] <- torch::torch_tensor(t(nuisance$mu[, , t])) # n x k
    }
    mu_a_star <- mu_star_nodes[[1]]$detach()$clone() # psi = ICE at node 1
    mu_a_star$requires_grad_(TRUE)
    beta_star <- beta

    for (tmle_iter in 1:tmle_maxiter) {
      max_eps <- 0

      for(t in tau:1) {
        mu_a_star <- mu_a_star$detach()$clone()
        mu_a_star$requires_grad_(TRUE)

        clever_K <- calculate_K(combined_loss, mu_a_star, Q_star, beta_star, design_matrix)

        cleverA_t <- calculate_clever_node(
          combined_loss, HA_node[, , t], mu_a_star, Q_star, beta_star, design_matrix
        )

        epsilon_t <- tmle_mle(
          p, node_fluctuation_model,
          mu_star_nodes[[t]],
          mu_star_nodes[[t + 1]],
          cleverA_t, clever_K, Q_star
        )

        max_eps <- max(max_eps, max(abs(as.numeric(epsilon_t))))
        cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {max_eps}\n\n"))

        # Update the node-t fit for each regime
        upd <- mu_star_nodes[[t]]$detach()$clone()
        for(j in 1:k) {
          if(tmle_linear == TRUE) {
            upd[, j] <- mu_star_nodes[[t]][, j] + cleverA_t[j, , ]$matmul(epsilon_t)
          }
          else {
            upd[, j] <- torch::torch_sigmoid(
              mu_star_nodes[[t]][, j]$logit() + cleverA_t[j, , ]$matmul(epsilon_t)
            )
          }
        }
        mu_star_nodes[[t]] <- upd$detach()

        # Fluctuate Q only once per outer iteration (e.g., at the final node t = 1)
        if(t == 1) {
          Q_star <- Q_fluctuation(epsilon_t, clever_K, Q_star)
        }
      }

      # After the backward sweep, node 1 holds the updated psi = ICE at L1
      mu_a_star <- mu_star_nodes[[1]]$detach()$clone()
      mu_a_star$requires_grad_(TRUE)
      beta_star <- B(combined_loss, mu_a_star$detach(), design_matrix, Q_star$detach())$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if (abs(max_eps) < 1e-2) break
    }

    # Final TMLE plug-in and EIF
    mu_a_star <- mu_star_nodes[[1]]$detach()$clone()
    mu_a_star$requires_grad_(TRUE)

    tmle_est <- B(combined_loss, mu_a_star$detach(), design_matrix, Q_star$detach())

    # Delta reuses the telescoping ICE residual, now computed from the updated (fluctuated) nodes
    resid_star <- array(0, dim = c(k, n, tau))
    for(t in 1:tau) {
      resid_star[, , t] <- t(as.matrix(mu_star_nodes[[t + 1]] - mu_star_nodes[[t]]))
    }
    resid_star <- resid_star / pi_cumprod * W
    Delta_star <- t(apply(resid_star, c(1, 2), sum))

    tmle_eif <- eif(
      combined_loss, mu_a_star, tmle_est, design_matrix, Q_star$detach(), torch::torch_tensor(Delta_star)
    )

    tmle_se <- apply(tmle_eif, 2, stats::sd) / sqrt(n)
    tmle_lower <- tmle_est + stats::qnorm(0.025) * tmle_se
    tmle_upper <- tmle_est + stats::qnorm(0.975) * tmle_se
  }

  res <- list(
    estimand = "longitudinal_dose_response",
    p = p,
    n = n,
    tau = tau,
    formula = formula,
    working_model = working_model,
    loss = loss,
    terms = terms,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    nuisance = nuisance,
    regimes = regimes,
    plugin = list(
      est = as.numeric(plugin)
    ),
    onestep = list(
      est = as.numeric(onestep_est$est),
      se = as.numeric(onestep_est$se),
      lower = as.numeric(onestep_est$lower),
      upper = as.numeric(onestep_est$upper),
      eif = onestep_est$eif,
      psi = as.numeric(psi),
      joint_draws = onestep_joint
    )
  )

  if (tmle == TRUE) {
    res$tmle <- list(
      est = as.numeric(tmle_est),
      se = as.numeric(tmle_se),
      lower = as.numeric(tmle_lower),
      upper = as.numeric(tmle_upper),
      eif = tmle_eif
    )
  }

  class(res) <- "automsm"

  res
}

#' Helper for retrieving all columns up to just before $A_t$
#' @noRd
history_cols <- function(t, Ls, As) {
  if (t == 1) {
    Ls[[1]]
  } else {
    unlist(c(Ls[1:t], As[1:(t - 1)]))
  }
}

#' Cumulative on-protocol indicator per regime
#'
#' Returns a k x n x tau array where entry \[j, i, t\] is 1 if subject i's
#' observed treatment history matches regime j through time t, and 0
#' once the subject first deviates from regime j (and thereafter).
#'
#' @param regimes k x tau matrix of regime trajectories in {0, 1}
#' @param A_obs n x tau matrix of observed treatments (0/1)
#' @noRd
on_protocol_weights <- function(regimes, A_obs) {
  tau <- ncol(A_obs)
  n <- nrow(A_obs)
  k <- nrow(regimes)
  stopifnot(ncol(regimes) == tau)

  W <- array(0, dim = c(k, n, tau))
  for(j in 1:k) {
    regime <- as.numeric(regimes[j, ])
    match_t <- sweep(A_obs, 2, regime, FUN = "==")
    on_protocol <- t(apply(match_t, 1, function(x) as.logical(cumprod(x))))
    W[j,,] <- on_protocol
  }

  W
}

#' Form cumulative propensity scores
#'
#' Returns a k x n x tau array where entry \[j, i, t\] is the cumulative
#' product prod_{s=1}^t g_s(a_s^{(j)} | H_s) of the per-time-point
#' propensity scores evaluated under regime j's trajectory.
#' @noRd
cumulative_propensity_scores <- function(regimes, pi1) {
  tau <- ncol(pi1)
  n <- nrow(pi1)
  k <- nrow(regimes)
  stopifnot(ncol(regimes) == tau)
  pi0 <- 1 - pi1
  pi_cumprod <- array(dim = c(k, n, tau))
  for (j in 1:nrow(regimes)) {
    regime <- regimes[j, , drop = FALSE]
    trt_mask <- matrix(regime == 1, nrow = n, ncol = tau, byrow = TRUE)
    ctl_mask <- matrix(regime == 0, nrow = n, ncol = tau, byrow = TRUE)
    pi_regime <- pi1 * trt_mask + pi0 * ctl_mask
    pi_cumprod[j, , ] <- t(apply(pi_regime, 1, cumprod))
  }

  pi_cumprod
}


#' @importFrom SuperLearner SuperLearner predict.SuperLearner
#' @noRd
estimate_longitudinal_dose_response_propensity_scores <- function(
  data,
  Ls,
  As,
  cv,
  folds,
  learners_trt,
  cv_control,
  epsilon
) {
  n <- nrow(data)
  tau <- length(As)

  pi_hat <- matrix(nrow = n, ncol = tau)

  # Propensity score fits
  for (fold in seq_along(cv)) {
    train <- cv[[fold]]$training_set
    valid <- cv[[fold]]$validation_set

    for (t in 1:tau) {
      Ht <- history_cols(t, Ls, As)
      pi_model <- SuperLearner::SuperLearner(
        Y = data[[As[t]]][train],
        X = data[train, Ht, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      pi_hat[valid, t] <- SuperLearner::predict.SuperLearner(
        pi_model,
        newdata = data[valid, c(Ht), drop = FALSE],
        onlySL = TRUE
      )$pred
      pi_hat[valid, t] <- bound(pi_hat[valid, t], 0, 1, epsilon)
    }
  }

  return(pi_hat)
}

#' @importFrom SuperLearner SuperLearner predict.SuperLearner
#' @noRd
estimate_longitudinal_dose_response_regressions <- function(
  data,
  Ls,
  As,
  Y,
  regimes,
  learners_outcome,
  cv,
  folds,
  outcome_type,
  cv_control,
  epsilon
) {
  n <- nrow(data)
  tau <- length(As)
  k <- nrow(regimes)

  # Recursion
  regress_and_predict <- function(t, regimes, train, valid) {
    # Base case: if t == tau, then pseudo-outcome is Y
    if (t == tau + 1) {
      return(list(list(
        regime = c(),
        train = matrix(data[[Y]][train], ncol = 1),
        valid = matrix(data[[Y]][valid], ncol = 1)
      )))
    }

    # We only need pseudo-outcomes for *unique* future treatment trajectories
    outcomes <- regress_and_predict(
      t + 1,
      unique(regimes[, (t + 1):tau, drop = FALSE]),
      train,
      valid
    )
    Ht_At <- c(history_cols(t, Ls, As), As[t])

    current_outcomes <- list()

    # Now fit new models to each of the unique pseudo-outcome sets
    for (po_index in seq_along(outcomes)) {
      outcome <- outcomes[[po_index]]

      outcome_family <- stats::gaussian()
      if (outcome_type == "binomial") {
        outcome_family <- stats::binomial
      }

      mu_model <- SuperLearner::SuperLearner(
        Y = outcome$train[, 1],
        X = data[train, Ht_At, drop = FALSE],
        SL.library = learners_outcome,
        family = if (t == tau) outcome_family else gaussian(),
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      for (a_t in unique(regimes[, 1, drop = FALSE][[As[t]]])) {
        newdata_train <- data[train, Ht_At, drop = FALSE]
        newdata_valid <- data[valid, Ht_At, drop = FALSE]
        newdata_train[[As[t]]] <- a_t
        newdata_valid[[As[t]]] <- a_t

        current_outcomes[[length(current_outcomes) + 1]] <- list(
          regime = c(a_t, outcome$regime),
          train = cbind(
            SuperLearner::predict.SuperLearner(
              mu_model,
              newdata = newdata_train,
              onlySL = TRUE
            )$pred,
            outcome$train
          ),
          valid = cbind(
            SuperLearner::predict.SuperLearner(
              mu_model,
              newdata = newdata_valid,
              onlySL = TRUE
            )$pred,
            outcome$valid
          )
        )
      }
    }

    return(current_outcomes)
  }

  results <- list()
  mu_valid <- array(dim = c(k, n, tau + 1))
  for (fold in seq_len(folds)) {
    train <- cv[[fold]]$training_set
    valid <- cv[[fold]]$validation_set
    results <- regress_and_predict(1, regimes, train, valid)

    for (j in seq_along(results)) {
      mu_valid[j, valid, ] <- results[[j]]$valid
    }
  }

  mu_valid
}


#' @importFrom  origami  make_folds
#' @importFrom  SuperLearner SuperLearner.CV.control predict.SuperLearner SuperLearner
#' @importFrom stats gaussian binomial
#' @noRd
estimate_longitudinal_dose_response_nuisance <- function(
  data,
  Ls,
  As,
  Y,
  regimes,
  learners_trt,
  learners_outcome,
  cv,
  outer_folds,
  inner_folds,
  outcome_type,
  epsilon = 1e-5
) {
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  pi_hat <- estimate_longitudinal_dose_response_propensity_scores(
    data,
    Ls,
    As,
    cv,
    outer_folds,
    learners_trt,
    cv_control,
    epsilon
  )

  mu_hat <- estimate_longitudinal_dose_response_regressions(
    data,
    Ls,
    As,
    Y,
    regimes,
    learners_outcome,
    cv,
    outer_folds,
    outcome_type,
    cv_control,
    epsilon
  )

  list(
    pi = pi_hat,
    mu = mu_hat
  )
}
