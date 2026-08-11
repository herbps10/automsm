#' Non-parametric marginal structural model estimation for dose response curves.
#'
#' @param data A \code{data.frame} containing the columns referenced by \code{X},
#'   \code{A}, and \code{Y}.
#' @param X A character vector giving the covariate column names in \code{data}.
#' @param A A string naming the treatment column in \code{data}.
#' @param Y A string naming the outcome column in \code{data}.
#' @param formula A model formula specifying the marginal structural working-model design matrix.
#' @param outcome_type Outcome type, either \code{"continuous"} or \code{"binomial"}.
#' @param loss The marginal structural model loss function measuring the fidelity of the
#'   working model to the target functional (e.g. \code{loss_squared_error}).
#' @param working_model The marginal structural model working model (e.g. \code{working_model_linear}).
#' @param learners_trt A character vector of \pkg{SuperLearner} libraries for estimating
#'   the propensity scores.
#' @param learners_outcome A character vector of \pkg{SuperLearner} libraries for
#'   estimating the outcome regressions.
#' @param outer_folds Number of folds in the outer cross-fitting loop.
#' @param inner_folds Number of folds for the inner \pkg{SuperLearner} cross-validation
#'   within each outer cross-fitting fold.
#' @param tmle Logical; whether to run the TMLE estimator.
#' @param tmle_maxiter Maximum number of TMLE iterations.
#' @param tmle_linear Logical; whether to use a linear TMLE fluctuation model
#'   (\code{TRUE}) or a logistic fluctuation model (\code{FALSE}).
#' @param bayes whether to run Bayesian TMLE estimator
#' @param bayes_draws number of MCMC samples
#' @param bayes_chains number of independent MCMC chains
#' @param bayes_prior prior to apply to marginal structural model parameters
#' @param epsilon Adjustment bounding estimated propensities/means away from 0 and 1.
#' @param nuisance Optional list of pre-computed nuisance parameters. If \code{NULL}
#'   (the default), nuisance parameters are estimated internally via cross-fitting.
#'
#' @return An object of class \code{"automsm"}: a list with components
#'   \describe{
#'     \item{estimand}{Character string, \code{"dose_response"}.}
#'     \item{p}{Number of working-model coefficients.}
#'     \item{n}{Sample size.}
#'     \item{formula}{The working-model formula used.}
#'     \item{working_model, loss}{The working model and loss function used.}
#'     \item{terms}{Character vector of working-model design-matrix term names.}
#'     \item{learners_trt, learners_outcome}{The \pkg{SuperLearner} libraries used.}
#'     \item{nuisance}{The (estimated or supplied) nuisance parameters.}
#'     \item{plugin}{A list with the plug-in piont estimate (\code{est}).}
#'     \item{onestep}{A list with the one-step point estimate (\code{est}), standard
#'       errors (\code{se}), confidence-interval bounds (\code{lower}, \code{upper}),
#'       the estimated efficient influence function (\code{eif}), the projected
#'       conditional means (\code{psi}), and joint draws (\code{joint_draws}).}
#'   }
#'
#' @seealso \code{\link{dose_response}} for the analogous estimator with a
#'  high-dimensional (categorical) point treatment.
#'
#' @importFrom stats sd model.matrix var qnorm dnorm
#' @importFrom torch torch_tensor torch_zeros torch_reshape nn_bce_with_logits_loss nn_mse_loss
#' @importFrom adaptMCMC MCMC
#'
#' @export
dose_response <- function(
  data,
  X,
  A,
  Y,
  formula,
  outcome_type = "binomial",
  loss = loss_squared_error,
  working_model = working_model_linear,
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm",
  outer_folds = 5,
  inner_folds = 5,
  tmle = TRUE,
  tmle_maxiter = 25,
  tmle_linear = TRUE,
  bayes = FALSE,
  bayes_draws = 1e3,
  bayes_chains = 4,
  bayes_prior = function(beta) {
    sum(stats::dnorm(as.numeric(beta), mean = 0, sd = 1, log = TRUE))
  },
  epsilon = 1e-5,
  nuisance = NULL
) {
  # ----- Argument checks -----
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)

  # Check A
  checkmate::assert_string(A)
  checkmate::assert_subset(A, choices = names(data))

  # Check X
  checkmate::assert_character(X, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(X, choices = names(data))

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
  if (is.null(nuisance)) {
    nuisance <- estimate_dose_response_nuisance(
      data,
      X,
      A,
      Y,
      learners_trt,
      learners_outcome,
      outer_folds,
      inner_folds,
      outcome_type,
      estimate_conditional_variance = bayes,
      epsilon = epsilon
    )
  }

  Xt <- torch::torch_tensor(as.matrix(data[, X, drop = FALSE]))
  Yt <- torch::torch_tensor(data[[Y]])

  As <- sort(unique(data[[A]]))
  K <- length(As)

  mat <- stats::model.matrix(formula, data = data)
  terms <- colnames(mat)
  p <- ncol(mat)

  design_matrix <- torch::torch_zeros(K * n * ncol(mat))$reshape(c(
    n,
    K,
    ncol(mat)
  ))

  for (k in 1:K) {
    datak <- data
    datak[[A]] <- As[k]
    mat <- model.matrix(formula, data = datak)
    terms <- colnames(mat)
    design_matrix[, k, ] <- mat
  }

  Delta <- function(H, Y, mu) {
    n <- nrow(H)
    H$mul((Y - mu)$reshape(c(n, 1)))
  }

  ##### Plugin estimator
  Q <- torch::torch_tensor(rep(1 / n, n))
  psi <- torch::torch_tensor(nuisance$mu_a, requires_grad = TRUE)
  plugin <- B(Lm(loss, working_model), psi, design_matrix, Q)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  H <- torch_zeros(c(n, K))
  HA <- torch_zeros(c(n, K))
  for (k in 1:K) {
    H[, k] <- (data[[A]] == As[k]) / nuisance$pi
    HA[, k] <- 1 / nuisance$pi_a[, k]
  }

  ##### One-step estimator
  onestep_est <- onestep(
    Lm(loss, working_model),
    psi,
    beta,
    design_matrix,
    Q,
    Delta(H, Yt, torch::torch_tensor(nuisance$mu))
  )

  # ----- TMLE -----
  if (tmle == TRUE) {
    # Calculate clever covariates for fluctuation model
    calculate_clever <- function(Lm, H, HA, psi, Q, beta, design_matrix, Minv) {
      n <- dim(design_matrix)[1]
      k <- dim(design_matrix)[2]
      p <- dim(design_matrix)[3]

      clever <- torch::torch_tensor(matrix(0, n, p))
      cleverA <- torch::torch_tensor(array(0, dim = c(n, k, p)))
      for (i in 1:n) {
        d <- Minv$matmul(grad_dL(
          Lm,
          psi[i, drop = FALSE],
          beta,
          design_matrix[i, , drop = FALSE]
        ))
        clever[i, ] <- torch::torch_reshape(d$matmul(H[i, ]), p)
        for (j in 1:k) {
          x <- torch::torch_zeros(k)
          x[k] <- HA[i, j]
          cleverA[i, k, ] <- torch::torch_reshape(d$matmul(x), p)
        }
      }
      list(
        clever = clever,
        cleverA = cleverA
      )
    }

    # TMLE fluctuation model
    if (tmle_linear == TRUE) {
      tmle_loss <- nn_mse_loss(reduction = "sum")
    } else {
      tmle_loss <- nn_bce_with_logits_loss(reduction = "sum")
    }

    tmle_fluctuation_model <- function(
      epsilon,
      mu,
      mu_a,
      clever,
      clever_K,
      Q,
      Y,
      design_matrix,
      condvar = NULL,
      bayes = FALSE
    ) {
      if (tmle_linear == TRUE) {
        mu <- mu + clever$clever$matmul(epsilon)
      } else {
        mu <- mu$logit() + clever$clever$matmul(epsilon)
      }

      Q <- Q_fluctuation(epsilon, clever_K, Q)

      if (bayes == FALSE) {
        target <- tmle_loss(mu, Y)
        target <- target - log(Q)$sum()
        target
      } else {
        if (tmle_linear == TRUE) {
          mu_a <- mu_a + clever$cleverA$matmul(epsilon)
        } else {
          mu_a <- torch::torch_sigmoid(mu_a$logit() + clever$cleverA$matmul(epsilon))
        }

        beta <- B(Lm(loss, working_model), mu_a, design_matrix, Q)

        if (tmle_linear == TRUE) {
          target <- as.numeric(
            -((mu - Y)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()
          ) +
            as.numeric(log(Q)$sum())
        } else {
          target <- -target
        }

        # Prior
        target <- target + bayes_prior(as.numeric(beta))
        jacobian <- torch::torch_zeros(c(p, p))
        J <- dB_dpsi(Lm(loss, working_model), mu_a, Q, design_matrix, beta)
        for (k in 1:K) {
          if (tmle_linear == TRUE) {
            jacobian <- jacobian +
              J[k, , ]$transpose(1, 2)$matmul(clever$cleverA[, k, ])
          } else {
            jacobian <- jacobian +
              J[k, , ]$transpose(1, 2)$matmul(
                (clever$cleverA[, k, ] * mu_a[, k]) * (1 - mu_a[, k, ])
              )
          }
        }

        jacobian <- jacobian +
          torch::torch_transpose(
            dB_dQ(Lm(loss, working_model), psi, Q, design_matrix, beta),
            1,
            2
          )$matmul(dQ_fluctuation_depsilon(epsilon, clever_K, Q))
        target <- target + log(abs(jacobian$det()))

        ret <- list(log.density = as.numeric(target), beta = as.numeric(beta))

        return(ret)
      }
    }

    Q_star <- Q
    mu_star <- torch::torch_tensor(nuisance$mu)
    mu_a_star <- torch::torch_tensor(nuisance$mu_a)
    beta_star <- beta
    epsilon_star <- rep(0, p)

    for (tmle_iter in 1:tmle_maxiter) {
      mu_a_star <- mu_a_star$detach()$clone()
      mu_a_star$requires_grad_(TRUE)

      Minv <- normalizing_matrix(Lm(loss, working_model), mu_a_star, beta_star, design_matrix, Q_star)

      clever_K <- calculate_K(
        Lm(loss, working_model),
        mu_a_star,
        Q_star,
        beta_star,
        design_matrix,
        Minv
      )

      clever <- calculate_clever(
        Lm(loss, working_model),
        H,
        HA,
        mu_a_star,
        Q_star,
        beta_star,
        design_matrix,
        Minv
      )

      epsilon_star <- tmle_mle(
        p,
        tmle_fluctuation_model,
        mu_star,
        mu_a_star,
        clever,
        clever_K,
        Q_star,
        Yt
      )

      m <- max(as.numeric(epsilon_star))
      #cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      if (tmle_linear == TRUE) {
        mu_star <- mu_star + clever$clever$matmul(epsilon_star)
      } else {
        mu_star <- torch::torch_sigmoid(
          mu_star$logit() + clever$clever$matmul(epsilon_star)
        )
      }

      mu_a_star <- mu_a_star$detach()$clone()
      for (k in 1:K) {
        if (tmle_linear == TRUE) {
          mu_a_star[, k] <- mu_a_star[, k] +
            clever$cleverA[, k, ]$matmul(epsilon_star)
        } else {
          mu_a_star[, k] <- torch::torch_sigmoid(
            mu_a_star[, k]$logit() + clever$cleverA[, k, ]$matmul(epsilon_star)
          )
        }
      }

      Q_star <- Q_fluctuation(epsilon_star, clever_K, Q_star)

      beta_star <- B(
        Lm(loss, working_model),
        mu_a_star$detach(),
        design_matrix,
        Q_star$detach()
      )$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if (abs(m) < 1e-2) {
        break
      }
    }

    mu_a_star <- mu_a_star$detach()$clone()
    mu_a_star$requires_grad_(TRUE)

    tmle_est <- B(
      Lm(loss, working_model),
      mu_a_star$detach(),
      design_matrix,
      Q_star$detach()
    )

    tmle_eif <- eif(
      Lm(loss, working_model),
      mu_a_star,
      tmle_est,
      design_matrix,
      Q_star$detach(),
      Delta(H, Yt, mu_star$detach())
    )
    tmle_se <- apply(tmle_eif, 2, stats::sd) / sqrt(n)
    tmle_lower <- tmle_est + stats::qnorm(0.025) * tmle_se
    tmle_upper <- tmle_est + stats::qnorm(0.975) * tmle_se

    tmle_beta_samples <- NULL
    tmle_acc_rate <- NULL

    if (bayes == TRUE) {
      tmle_beta_samples <- array(dim = c(bayes_chains, bayes_draws, p))
      tmle_acc_rate <- 0

      for (chain in 1:bayes_chains) {
        log_dens <- function(epsilon) {
          tmle_fluctuation_model(
            torch::torch_tensor(epsilon),
            mu_star,
            mu_a_star,
            clever,
            clever_K,
            Q_star,
            Yt,
            design_matrix,
            condvar = torch::torch_tensor(nuisance$condvar),
            bayes = TRUE
          )
        }

        mcmc <- adaptMCMC::MCMC(
          log_dens,
          n = bayes_draws,
          init = as.numeric(epsilon_star),
          adapt = TRUE,
          acc.rate = 0.3,
          scale = rep(1e-3, p)
        )

        tmle_beta_samples[chain, , ] <- matrix(
          unlist(mcmc$extra.values),
          ncol = p,
          nrow = bayes_draws,
          byrow = TRUE
        )
        tmle_acc_rate <- tmle_acc_rate + 1 / bayes_chains * mcmc$acceptance.rate
      }
    }
  }

  res <- list(
    estimand = "dose_response",
    p = p,
    n = n,
    formula = formula,
    working_model = working_model,
    loss = loss,
    terms = terms,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    nuisance = nuisance,
    plugin = list(
      est = as.numeric(plugin)
    ),
    onestep = list(
      est = as.numeric(onestep_est$est),
      se = as.numeric(onestep_est$se),
      lower = as.numeric(onestep_est$lower),
      upper = as.numeric(onestep_est$upper),
      eif = onestep_est$eif,
      psi = as.numeric(psi)
    )
  )

  if (tmle == TRUE) {
    res$tmle <- list(
      est = as.numeric(tmle_est),
      se = as.numeric(tmle_se),
      lower = as.numeric(tmle_lower),
      upper = as.numeric(tmle_upper),
      eif = tmle_eif,
      samples = tmle_beta_samples,
      acc_rate = tmle_acc_rate
    )
  }

  class(res) <- "automsm"

  res
}

#' @importFrom  origami  make_folds
#' @importFrom  SuperLearner SuperLearner.CV.control predict.SuperLearner SuperLearner
#' @importFrom stats gaussian binomial
#' @noRd
estimate_dose_response_nuisance <- function(
  data,
  X,
  A,
  Y,
  learners_trt,
  learners_outcome,
  outer_folds,
  inner_folds,
  outcome_type,
  estimate_conditional_variance = FALSE,
  epsilon = 1e-5
) {
  n <- nrow(data)
  As <- sort(unique(data[[A]]))
  k <- length(As)

  pi_a_hat <- mu_a_hat <- matrix(0, n, k)
  pi_hat <- mu_hat <- condvar_hat <- numeric(n)

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  if (outer_folds == 1) {
    cv[[1]]$training_set <- cv[[1]]$validation_set
  }
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  outcome_family <- stats::gaussian()
  if (outcome_type == "binomial") {
    outcome_family <- stats::binomial()
  }

  for (fold in seq_along(cv)) {
    training <- cv[[fold]]$training_set
    validation <- cv[[fold]]$validation_set

    for (a_index in seq_along(As)) {
      pi_model <- SuperLearner::SuperLearner(
        Y = as.numeric(data[[A]][training] == As[a_index]),
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      pi_a_hat[validation, a_index] <- SuperLearner::predict.SuperLearner(
        pi_model,
        newdata = data[validation, c(X), drop = FALSE],
        onlySL = TRUE
      )$pred
      pi_a_hat[validation, a_index] <- bound(
        pi_a_hat[validation, a_index],
        0,
        1,
        epsilon
      )
    }

    mu_model <- SuperLearner::SuperLearner(
      Y = data[[Y]][training],
      X = data[training, c(X, A), drop = FALSE],
      SL.library = learners_outcome,
      family = outcome_family,
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    for (a_index in seq_along(As)) {
      newdata <- data[validation, c(X, A)]
      newdata[[A]] <- As[a_index]
      mu_a_hat[validation, a_index] <- SuperLearner::predict.SuperLearner(
        mu_model,
        newdata = newdata,
        onlySL = TRUE
      )$pred
      if (outcome_type == "binomial") {
        mu_a_hat[validation, a_index] <- bound(
          mu_a_hat[validation, a_index],
          0,
          1,
          epsilon
        )
      }
    }

    if (estimate_conditional_variance == TRUE) {
      yhat <- SuperLearner::predict.SuperLearner(
        mu_model,
        newdata = data[training, c(X, A)],
        onlySL = TRUE
      )$pred
      y2 <- (data[[Y]][training] - yhat)^2

      yvar_model <- SuperLearner::SuperLearner(
        Y = y2,
        X = data[training, c(X, A), drop = FALSE],
        SL.library = learners_outcome,
        family = outcome_family,
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      condvar_hat[validation] <- SuperLearner::predict.SuperLearner(
        yvar_model,
        newdata = data[validation, c(X, A)],
        onlySL = TRUE
      )$pred
      condvar_hat[validation] <- ifelse(
        condvar_hat[validation] < epsilon,
        epsilon,
        condvar_hat[validation]
      )
    }
  }

  for (a_index in seq_along(As)) {
    ind <- data[[A]] == As[a_index]
    mu_hat[ind] <- mu_a_hat[ind, a_index]
    pi_hat[ind] <- pi_a_hat[ind, a_index]
  }

  list(
    pi_a = pi_a_hat,
    mu_a = mu_a_hat,
    mu = mu_hat,
    pi = pi_hat,
    condvar = condvar_hat
  )
}

