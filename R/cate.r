#' Estimate Marginal Structural Model for the Conditional Average Treatment Effect (CATE)
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
#' #' @return An object of class \code{"automsm"}: a list with components
#'   \describe{
#'     \item{estimand}{Character string, \code{"cate"}.}
#'     \item{p}{Number of working-model coefficients.}
#'     \item{n}{Sample size.}
#'     \item{formula}{The working-model formula used.}
#'     \item{working_model, loss}{The working model and loss function used.}
#'     \item{terms}{Character vector of working-model design-matrix term names.}
#'     \item{learners_trt, learners_outcome}{The \pkg{SuperLearner} libraries used.}
#'     \item{nuisance}{The (estimated or supplied) nuisance parameters.}
#'     \item{plugin}{A list with the plug-in point estimate (\code{est}).}
#'     \item{onestep}{A list with the one-step point estimate (\code{est}), standard
#'       errors (\code{se}), confidence-interval bounds (\code{lower}, \code{upper}),
#'       the estimated efficient influence function (\code{eif}), and joint draws
#'       (\code{joint_draws}).}
#'   }

#' @importFrom torch torch_tensor torch_zeros torch_reshape nn_bce_with_logits_loss nn_mse_loss
#' @importFrom stats dnorm qnorm
#' @importFrom adaptMCMC MCMC
#'
#' @export
cate <- function(
  data,
  X,
  A,
  Y,
  formula,
  outcome_type = "binomial",
  loss = loss_weighted_sum(loss_squared_error),
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
  validate_msm_arguments(
    data, X, A, Y, formula,
    outcome_type, loss, working_model,
    learners_trt, learners_outcome,
    outer_folds, inner_folds,
    tmle, tmle_maxiter, tmle_linear,
    bayes, bayes_draws, bayes_chains, bayes_prior,
    epsilon, nuisance
  )

  n <- nrow(data)

  if (is.null(nuisance)) {
    nuisance <- estimate_cate_nuisance(
      data,
      X,
      A,
      Y,
      learners_trt,
      learners_outcome,
      outer_folds,
      inner_folds,
      outcome_type,
      estimate_conditional_variance = bayes && tmle_linear == TRUE,
      epsilon = epsilon
    )
  }

  psi <- torch::torch_tensor(nuisance$mu1 - nuisance$mu0, requires_grad = TRUE)$reshape(c(n, 1))

  # Design matrix
  mat <- model.matrix(formula, data = data)
  terms <- colnames(mat)
  design_matrix <- torch::torch_tensor(mat)$reshape(c(n, 1, ncol(mat)))
  p <- ncol(mat)

  Yt <- torch::torch_tensor(data[[Y]], dtype = torch::torch_float())
  H0 <- torch_tensor(-1 / (1 - nuisance$pi))
  H1 <- torch_tensor(1 / nuisance$pi)
  H <- torch_tensor(ifelse(data[[A]] == 1, 1 / nuisance$pi, -1 / (1 - nuisance$pi)))
  Q <- torch::torch_tensor(rep(1 / n, n))

  Delta <- (H * (Yt - torch::torch_tensor(nuisance$mu)))$reshape(c(n, 1))

  base <- estimate_plugin_and_onestep(
    loss,
    working_model,
    design_matrix,
    psi,
    Q,
    Delta
  )

  ##### TMLE

  tmle_res <- NULL
  bayes_tmle_res <- NULL

  if (tmle == TRUE) {
    # Calculate clever covariates for fluctuation model
    calculate_clever <- function(Lm, H, H0, H1, psi, Q, beta, design_matrix, Minv) {
      n <- dim(design_matrix)[1]
      p <- dim(design_matrix)[3]

      clever <- torch::torch_tensor(matrix(0, n, p))
      clever0 <- torch::torch_tensor(matrix(0, n, p))
      clever1 <- torch::torch_tensor(matrix(0, n, p))

      for (i in 1:n) {
        d <- Minv$matmul(grad_dL(Lm, psi[i, drop = FALSE], beta, design_matrix[i, , drop = FALSE]))
        clever[i, ] <- torch::torch_reshape(H[i] * d, p)
        clever0[i, ] <- torch::torch_reshape(H0[i] * d, p)
        clever1[i, ] <- torch::torch_reshape(H1[i] * d, p)
      }
      list(
        clever = clever,
        clever0 = clever0,
        clever1 = clever1
      )
    }

    # TMLE fluctuation model
    if (tmle_linear == TRUE) {
      tmle_loss <- torch::nn_mse_loss(reduction = "sum")
    } else {
      tmle_loss <- torch::nn_bce_with_logits_loss(reduction = "sum")
    }

    tmle_fluctuation_model <- function(
      epsilon,
      mu,
      mu0,
      mu1,
      clever,
      clever0,
      clever1,
      K,
      Q,
      Y,
      design_matrix,
      Lm = NULL,
      condvar = NA,
      bayes = FALSE
    ) {
      # Fluctuate covariate distribution
      #Q <- Q_fluctuation(epsilon, K, Q)

      # Fluctuate outcome regression
      if (tmle_linear == TRUE) {
        mu <- mu[, 1] + clever$matmul(epsilon)
        mu0 <- mu0[, 1] + clever0$matmul(epsilon)
        mu1 <- mu1[, 1] + clever1$matmul(epsilon)
        target <- tmle_loss(mu, Y)
      } else {
        mu_logit <- mu[, 1]$logit() + clever$matmul(epsilon)
        mu0 <- torch::torch_sigmoid(mu0[,1]$logit() + clever0$matmul(epsilon))
        mu1 <- torch::torch_sigmoid(mu1[,1]$logit() + clever1$matmul(epsilon))
        target <- tmle_loss(mu_logit[, 1], Y)
      }

      # Combined loss function
      target <- target - Q$log()$sum()

      if (bayes == FALSE) {
        return(target)
      } else {
        psi <- (mu1 - mu0)$reshape(c(n, 1))
        beta <- B(
          Lm(loss, working_model),
          psi$detach()$clone(),
          design_matrix,
          Q$detach()$clone()
        )

        if (tmle_linear == TRUE) {
          target <- as.numeric(
            -((mu - Y)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()
          ) +
            as.numeric(log(Q)$sum())
        } else {
          target <- -as.numeric(target)
        }

        # Prior
        target <- target + bayes_prior(as.numeric(beta))
        # Jacobian adjustment
        jacobian <- dB_dpsi(Lm(loss, working_model), psi, Q, design_matrix, beta)$transpose(2, 3)
        if (tmle_linear == TRUE) {
          jacobian <- jacobian$matmul(clever1 - clever0)
        } else {
          jacobian <- jacobian$matmul(
            clever1 *
              mu1$reshape(c(n, 1)) *
              (1 - mu1$reshape(c(n, 1))) -
              clever0 * mu0$reshape(c(n, 1)) * (1 - mu0$reshape(c(n, 1)))
          )
        }

        jacobian <- jacobian +
          torch::torch_transpose(
            dB_dQ(Lm(loss, working_model), psi, Q, design_matrix, beta),
            1,
            2
          )$matmul(dQ_fluctuation_depsilon(epsilon, K, Q))
        target <- target + log(abs(jacobian$det()))

        ret <- list(log.density = as.numeric(target), beta = as.numeric(beta))

        return(ret)
      }
    }

    Q_star <- Q
    mu_star <- torch::torch_tensor(nuisance$mu)$reshape(c(n, 1))
    mu0_star <- torch::torch_tensor(nuisance$mu0)$reshape(c(n, 1))
    mu1_star <- torch::torch_tensor(nuisance$mu1)$reshape(c(n, 1))
    beta_star <- torch::torch_tensor(base$plugin$est, requires_grad = TRUE)
    epsilon_star <- rep(0, p)

    converged <- TRUE

    for (tmle_iter in 1:tmle_maxiter) {
      psi_star <- (mu1_star - mu0_star)$detach()$clone()
      psi_star$requires_grad_(TRUE)

      Minv <- normalizing_matrix(Lm(loss, working_model), psi_star, beta_star, design_matrix, Q_star)

      K <- calculate_K(
        Lm(loss, working_model),
        psi_star,
        beta_star,
        design_matrix,
        Minv
      )

      clever <- calculate_clever(
        Lm(loss, working_model),
        H,
        H0,
        H1,
        psi_star,
        Q_star,
        beta_star,
        design_matrix,
        Minv
      )

      epsilon_star <- tmle_mle(p, function(epsilon) {
        tmle_fluctuation_model(
          epsilon, mu_star, mu0_star, mu1_star,
          clever$clever, clever$clever0, clever$clever1,
          K, Q_star, Yt, design_matrix
        )
      })

      if (any(is.nan(as.numeric(epsilon_star)))) {
        warning("TMLE failed to converge.")
        converged <- FALSE
        break
      }

      m <- max(abs(as.numeric(epsilon_star)))
      #cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      if (tmle_linear == TRUE) {
        mu_star[, 1]  <- mu_star[, 1]  + clever$clever$matmul(epsilon_star)
        mu0_star[, 1] <- mu0_star[, 1] + clever$clever0$matmul(epsilon_star)
        mu1_star[, 1] <- mu1_star[, 1] + clever$clever1$matmul(epsilon_star)
      } else {
        mu_star[, 1]  <- torch::torch_sigmoid(mu_star[, 1]$logit()  + clever$clever$matmul(epsilon_star))
        mu0_star[, 1] <- torch::torch_sigmoid(mu0_star[, 1]$logit() + clever$clever0$matmul(epsilon_star))
        mu1_star[, 1] <- torch::torch_sigmoid(mu1_star[, 1]$logit() + clever$clever1$matmul(epsilon_star))
      }

      Q_star <- Q_fluctuation(epsilon_star, K, Q_star)

      beta_star <- B(
        Lm(loss, working_model),
        psi_star,
        design_matrix,
        Q_star
      )$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if (abs(m) < 1e-3) {
        break
      }
    }

    tmle_est <- tmle_se <- tmle_lower <- tmle_upper <- rep(NA, p)
    tmle_eif <- matrix(NA, ncol = p, nrow = n)
    if (converged == TRUE) {
      tmle_est <- B(Lm(loss, working_model), psi_star, design_matrix, Q_star)
      tmle_eif <- eif(
        Lm(loss, working_model),
        psi_star,
        tmle_est,
        design_matrix,
        Q_star,
        (H * (Yt - mu_star[, 1]))$reshape(c(n, 1))
      )
      tmle_se <- apply(tmle_eif, 2, sd) / sqrt(n)
      tmle_lower <- tmle_est + stats::qnorm(0.025) * tmle_se
      tmle_upper <- tmle_est + stats::qnorm(0.975) * tmle_se
    }

    tmle_res <- list(
      est = as.numeric(tmle_est),
      se = as.numeric(tmle_se),
      lower = as.numeric(tmle_lower),
      upper = as.numeric(tmle_upper),
      eif = tmle_eif
    )

    tmle_beta_samples <- NULL
    tmle_acc_rate <- NULL

    if (bayes == TRUE) {
      tmle_beta_samples <- array(dim = c(bayes_chains, bayes_draws, p))
      tmle_acc_rate <- 0

      for (chain in 1:bayes_chains) {
        cat("Chain: ", chain, "\n\n")

        log_dens <- function(epsilon) {
          tmle_fluctuation_model(
            torch::torch_tensor(epsilon),
            mu_star,
            mu0_star,
            mu1_star,
            clever$clever,
            clever$clever0,
            clever$clever1,
            K,
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

    bayes_tmle_res <- list(
      samples = tmle_beta_samples,
      acc_rate = tmle_acc_rate
    )
  }

  assemble_result(
    "cate",
    p,
    n,
    formula,
    working_model,
    loss,
    terms,
    learners_trt,
    learners_outcome,
    nuisance,
    plugin = list(est = as.numeric(base$plugin$est)),
    onestep = base$onestep,
    tmle_res,
    bayes_tmle_res
  )
}

#' @noRd
estimate_cate_nuisance <- function(
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
  nuis <- estimate_point_nuisance(
    data, X, A, Y,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    outcome_type = outcome_type,
    levels = c(0, 1),
    propensity = "binary",
    estimate_conditional_variance = estimate_conditional_variance,
    epsilon = epsilon
  )

  list(
    pi      = nuis$pi_a[, 2],
    mu0     = nuis$mu_a[, 1],
    mu1     = nuis$mu_a[, 2],
    mu      = nuis$mu,
    condvar = nuis$condvar
  )
}


