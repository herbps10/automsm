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
#' @param regimes The set of treatment trajectories \eqn{\bar{a} \in \{0, 1\}^\tau}
#'   over which the NP-MSM is defined. Supplied as a \code{data.frame} (or matrix)
#'   with \eqn{K} rows and \eqn{\tau} columns, where each row is one trajectory
#'   and the columns are named to match \code{As}.
#' @param summary_measures Optional function mapping the \code{regimes} \code{data.frame}
#'   to a \code{data.frame} of regime-level summary measures (one row per regime),
#'   used to build the marginal structural working-model design. For example,
#'   \code{function(regimes) data.frame(V = rowSums(regimes))} computes
#'   the cumulative treatment duration \eqn{v(\bar{a}) = \sum_t a_t}. When
#'   supplied, the design for each regime is constructed by
#'   evaluating \code{formula} on a data frame combining the observation-level
#'   covariates in \code{data} with that regime's summary measures. This
#'   permits formulas that mix observation-level effect modifiers and
#'   regime-level summaries. If \code{NULL} (the default), no summary
#'   measures are constructed. Note the choice of summary measures
#'   is part of the estimand definition, and different summaries define
#'   different target parameters (Petersen et al., 2014).
#' @param outcome_type Outcome type, either \code{"continuous"} or \code{"binomial"}.
#' @param loss The marginal structural model loss function measuring the fidelity of the
#'   working model to the target functional (e.g. \code{loss_weighted_sum}).
#' @param working_model The marginal structural model working model (e.g. \code{working_model_linear}).
#' @param p Integer giving the dimension of the working-model parameter
#'   \eqn{\beta} (the number of working-model coefficients). If \code{NULL}
#'   (the default), \code{p} is inferred from the number of columns in the
#'   design matrix produced by \code{formula}. When a working model whose
#'   coefficient dimension does not match the \code{formula}-based design, \code{p}
#'   should be supplied explicitly to match the working model's true number of coefficients.
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
  regimes,
  summary_measures = NULL,
  outcome_type = "binomial",
  loss = loss_squared_error,
  working_model = working_model_linear,
  p = NULL,
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
  validate_longitudinal_arguments(
    data, Ls, As, Y, formula, summary_measures,
    outcome_type, loss, working_model, p,
    learners_trt, learners_outcome,
    outer_folds, inner_folds,
    tmle, tmle_maxiter, tmle_linear,
    epsilon, nuisance
  )

  n <- nrow(data)
  tau <- length(As)
  k <- nrow(regimes)
  stopifnot(identical(colnames(regimes), As))

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

  if(!is.null(summary_measures)) {
    # ---- Convenience path: build k x n x p tensor from formula + summaries ----
    sm <- as.data.frame(summary_measures(regimes))
    checkmate::assert_data_frame(sm, nrows = k)
    if(any(colnames(sm) %in% names(data))) {
      stop("summary_measures returned column name(s) that collide with `data`: ", paste(intersect(colnames(sm), names(data)), collapse = ", "))
    }

    dm <- build_design_tensor(
      formula, data, K = k,
      mutate = function(d, j) {
        for(col in colnames(sm)) d[[col]] <- sm[j, col]
        d
      },
      p = p
    )
  }
  else {
    dm <- build_design_tensor(formula, data, K = k, mutate = NULL, p = p)
  }

  pi_cumprod <- cumulative_propensity_scores(regimes, nuisance$pi)
  W <- on_protocol_weights(regimes, as.matrix(data[, As, drop = FALSE]))
  HA_node <- W / pi_cumprod
  stopifnot(all(is.finite(HA_node)))

  problem <- new_msm_problem(
    estimand = "longitudinal_dose_response", K = k, d = dm$d, p = dm$p, tau = tau,
    design_matrix = dm$design_matrix,
    Q0 = torch::torch_tensor(rep(1 / n, n)),
    Yt = data[[Y]],
    Lm_fn = Lm(loss, working_model),
    loss = loss, working_model = working_model,
    formula = formula, terms = dm$terms,
    outcome_type = outcome_type, nuisance = nuisance,
    aux = list(
      HA_node = torch::torch_tensor(HA_node),
      pi_cumprod = pi_cumprod,
      W = W,
      regimes = regimes,
      Ls = Ls,
      As = As,
      cv = cv
    )
  )

  spec <- msm_spec_longitudinal_dose_response(tmle_linear = tmle_linear)

  res <- fit_msm(
    problem, spec,
    tmle = tmle_control(tmle, tmle_maxiter, tmle_linear),
    bayes = bayes_control(enabled = FALSE),
    onestep = onestep_control(1e3, 1)
  )

  assemble_result(
    "longitudinal_dose_response", problem$p, problem$d, problem$n, formula, working_model,
    loss, problem$terms, learners_trt, learners_outcome, nuisance,
    plugin = res$base$plugin, onestep = res$base$onestep,
    tmle = res$tmle,
    bayes_tmle = NULL,
    tau = tau, regimes = regimes
  )
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
  learners_trt,
  cv_control,
  epsilon
) {
  n <- nrow(data)
  tau <- length(As)
  pi_hat <- matrix(nrow = n, ncol = tau)
  for (fold in seq_along(cv)) {
    train <- cv[[fold]]$training_set
    valid <- cv[[fold]]$validation_set
    for (t in seq_len(tau)) {
      Ht <- history_cols(t, Ls, As)

      res <- sl_fit_predict(
        Y = data[[As[t]]][train],
        X = data[train, Ht, drop = FALSE],
        newdata = list(valid = data[valid, Ht, drop = FALSE]),
        SL.library = learners_trt,
        family = "binomial",
        cv_control = cv_control,
        bounds = c(0, 1),
        epsilon = epsilon
      )

      pi_hat[valid, t] <- res$pred$valid
    }
  }

  stopifnot(!anyNA(pi_hat))
  pi_hat
}

#' Estimate the sequential (ICE) outcome regressions for a longitudinal NP-MSM
#'
#' @details
#' The pseudo-outcomes at node \code{t} depends on the treatment trajectory only
#' through its suffix \eqn{(a_t, \dots, a_\tau)}. The recursion works over the
#' distinct suffixes actually present in \code{regimes}, matching each
#' to its parent suffix by key, so only the requested regimes are produced
#' (and slice \code{j} of the returned array corresponds to row \code{j} of
#' \code{regimes}). Each regression is fit pooling over all observed treatment
#' values and then evaluated at the regime's value of \eqn{A_t}, so rules
#' with little support still borrow strength from the full sample.
#'
#' @importFrom stats gaussian
#' @noRd
estimate_longitudinal_dose_response_regressions <- function(
  data,
  Ls,
  As,
  Y,
  regimes,
  learners_outcome,
  cv,
  outcome_type,
  outcome_family,
  cv_control,
  epsilon
) {
  n <- nrow(data)
  tau <- length(As)
  k <- nrow(regimes)

  bounds <- if(outcome_type == "binomial") c(0, 1) else NULL

  # Distinct suffixes requested at each level
  suffixes <- vector("list", tau + 1L)
  for(t in seq_len(tau)) {
    suffixes[[t]] <- unique(regimes[, As[t:tau], drop = FALSE])
  }
  suffixes[[tau + 1L]] <- regimes[1L, integer(0), drop = FALSE]

  # Recursion over suffixes
  regress_and_predict <- function(t, train, valid) {
    # Base case: if t == tau, then pseudo-outcome is Y
    if (t == tau + 1L) {
      return(list(
        keys = regime_key(suffixes[[tau + 1L]]),
        fits = list(list(
        train = matrix(data[[Y]][train], ncol = 1),
        valid = matrix(data[[Y]][valid], ncol = 1)
        ))
      ))
    }

    parent <- regress_and_predict(t + 1L, train, valid)

    needed <- suffixes[[t]]
    keys_t <- regime_key(needed)
    parent_idx <- match(regime_key(needed[, -1L, drop = FALSE]), parent$keys)
    stopifnot(!anyNA(parent_idx))

    a_t_vals <- needed[[As[t]]]
    Ht_At <- c(history_cols(t, Ls, As), As[t])
    Xtrain <- data[train, Ht_At, drop = FALSE]
    fam <- if(t == tau) outcome_family else stats::gaussian()

    out <- vector("list", nrow(needed))

    # One fit per distinct parent pseudo-outcome; reuse it across the a_t
    # values required by the regimes sharing that parent.
    for (pj in unique(parent_idx)) {
      po <- parent$fits[[pj]]

      mu_fit <- sl_fit(
        Y = po$train[, 1L],
        X = Xtrain,
        SL.library = learners_outcome,
        family = fam,
        cv_control = cv_control
      )

      for (j in which(parent_idx == pj)) {
        nd_train <- data[train, Ht_At, drop = FALSE]
        nd_valid <- data[valid, Ht_At, drop = FALSE]
        nd_train[[As[t]]] <- a_t_vals[j]
        nd_valid[[As[t]]] <- a_t_vals[j]

        out[[j]] <- list(
          train = cbind(sl_predict(mu_fit, nd_train, bounds, epsilon), po$train),
          valid = cbind(sl_predict(mu_fit, nd_valid, bounds, epsilon), po$valid)
        )
      }
    }

    list(keys = keys_t, fits = out)
  }

  regime_keys <- regime_key(regimes)
  mu_valid <- array(NA_real_, dim = c(k, n, tau + 1L))

  for (fold in seq_along(cv)) {
    train <- cv[[fold]]$training_set
    valid <- cv[[fold]]$validation_set

    res <- regress_and_predict(1L, train, valid)

    idx <- match(regime_keys, res$keys)
    stopifnot(!anyNA(idx))

    for (j in seq_len(k)) {
      mu_valid[j, valid, ] <- res$fits[[idx[j]]]$valid
    }
  }

  stopifnot(!anyNA(mu_valid))

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
  setup <- nuisance_setup(nrow(data), outer_folds, inner_folds, outcome_type, cv)

  pi_hat <- estimate_longitudinal_dose_response_propensity_scores(
    data,
    Ls,
    As,
    setup$cv,
    learners_trt,
    setup$cv_control,
    epsilon
  )

  mu_hat <- estimate_longitudinal_dose_response_regressions(
    data,
    Ls,
    As,
    Y,
    regimes,
    learners_outcome,
    setup$cv,
    outcome_type,
    setup$outcome_family,
    setup$cv_control,
    epsilon
  )

  list(
    pi = pi_hat,
    mu = mu_hat
  )
}

#' Stable key for a set of treatment trajectories (or suffixes thereof)
#'
#' Handles the zero-column case (the empty suffix at node \eqn{tau + 1})
#' @noRd
regime_key <- function(x) {
  x <- as.matrix(x)
  if(ncol(x) == 0L) return(rep("", nrow(x)))
  apply(x, 1L, paste, collapse = "\r")
}
