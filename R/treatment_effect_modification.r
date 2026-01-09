estimate_treatment_effect_modification_nuisance <- function(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type, estimate_conditional_variance = FALSE, epsilon = 1e-5) {
  n <- nrow(data)
  data0 <- data1 <- data
  data0[[A]] <- 0
  data1[[A]] <- 1
  pi_hat <- mu0_hat <- mu1_hat <- condvar_hat <- numeric(n)

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  if(outer_folds == 1) cv[[1]]$training_set <- cv[[1]]$validation_set
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  outcome_family <- stats::gaussian()
  if(outcome_type == "binomial") outcome_family <- stats::binomial()

  for(fold in seq_along(cv)) {
    training   <- cv[[fold]]$training_set
    validation <- cv[[fold]]$validation_set

    pi_model <- SuperLearner::SuperLearner(
      Y = data[[A]][training],
      X = data[training, X, drop = FALSE],
      SL.library = learners_trt,
      family = "binomial",
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    mu_model <- SuperLearner::SuperLearner(
      Y = data[[Y]][training],
      X = data[training, c(X, A), drop = FALSE],
      SL.library = learners_outcome,
      family = outcome_family,
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    pi_hat[validation]  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data[validation, c(X), drop = FALSE], onlySL = TRUE)$pred
    mu0_hat[validation] <- SuperLearner::predict.SuperLearner(mu_model, newdata = data0[validation, c(X, A)], onlySL = TRUE)$pred
    mu1_hat[validation] <- SuperLearner::predict.SuperLearner(mu_model, newdata = data1[validation, c(X, A)], onlySL = TRUE)$pred

    if(estimate_conditional_variance == TRUE) {
      yhat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data[training, c(X, A)], onlySL = TRUE)$pred
      y2 <- (data[[Y]][training] - yhat)^2

      yvar_model <- SuperLearner::SuperLearner(
        Y = y2,
        X = data[training, c(X, A), drop = FALSE],
        SL.library = learners_outcome,
        family = outcome_family,
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      condvar_hat[validation] <- SuperLearner::predict.SuperLearner(yvar_model, newdata = data[validation, c(X, A)], onlySL = TRUE)$pred
      condvar_hat[validation] <- ifelse(condvar_hat[validation] < 0, 0, condvar_hat[validation])
    }
  }

  if(outcome_type == "binomial") {
    mu0_hat[mu0_hat >= 1 - epsilon] <- 1 - epsilon
    mu0_hat[mu0_hat <= epsilon] <- epsilon

    mu1_hat[mu1_hat >= 1 - epsilon] <- 1 - epsilon
    mu1_hat[mu1_hat <= epsilon] <- epsilon
  }

  mu_hat <- ifelse(data[[A]] == 1, mu1_hat, mu0_hat)

  list(
    pi  = pi_hat,
    mu0 = mu0_hat,
    mu1 = mu1_hat,
    mu  = mu_hat,
    condvar = condvar_hat
  )
}

#' Estimate Marginal Structural Model for Treatment Effect Modifiers
#'
#' @param data data frame
#' @param X covariate columns
#' @param A treatment column
#' @param Y outcome column
#' @param formula marginal structural model design matrix formula
#' @param loss marginal structural model loss function
#' @param working model marginal structural model working model
#' @param learners_trt SuperLearner libraries for estimating propensity score
#' @param learners_outcome SuperLearner libraries for estimating outcome regression
#' @param outer_folds number of folds in outer cross-fitting loop
#' @param inner_folds number of folds for inner SuperLearner cross-validation within each outer cross-fitting loop
#' @param outcome_type outcome type ("continuous" or "binomial")
#' @param tmle whether to run TMLE estimator (TRUE/FALSE)
#' @param tmle_maxiter maximum number of TMLE iterations
#' @param tmle_linear whether to use linear TMLE fluctuation model (TRUE) or logistic fluctuation model (FALSE)
#' @param bayes whether to run Bayesian TMLE estimator
#' @param bayes_draws number of MCMC samples
#' @param bayes_prior prior to apply to marginal structural model parameters
#'
#' @export
treatment_effect_modification <- function(
    data,
    X,
    A,
    Y,
    formula,
    loss = loss_squared_error,
    working_model = working_model_linear,
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    outer_folds = 5,
    inner_folds = 5,
    outcome_type = "binomial",
    tmle = TRUE,
    tmle_maxiter = 25,
    tmle_linear = TRUE,
    bayes = FALSE,
    bayes_draws = 1e3,
    bayes_prior = \(beta) sum(dnorm(as.numeric(beta), mean = 0, sd = 1, log = TRUE)),
    epsilon = 1e-5,
    nuisance = NULL
) {
  n <- nrow(data)
  #data <- data[, c(X, A, Y)]
  if(is.null(nuisance)) {
    nuisance <- estimate_treatment_effect_modification_nuisance(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type, estimate_conditional_variance = bayes && tmle_linear == TRUE, epsilon = epsilon)
  }

  Xt <- torch::torch_tensor(as.matrix(data[, X, drop = FALSE]))
  Yt <- torch::torch_tensor(data[[Y]], dtype = torch::torch_float())

  mat <- model.matrix(formula, data = data)
  terms <- colnames(mat)
  design_matrix <- torch::torch_tensor(mat)$reshape(c(n, ncol(mat)))
  p <- ncol(design_matrix)

  Delta <- \(H, Y, mu) H * (Y - mu)

  ##### Plugin estimator
  Q <- torch::torch_tensor(rep(1 / n, n))
  psi <- torch::torch_tensor(nuisance$mu1 - nuisance$mu0, requires_grad = TRUE)
  plugin <- B(Lm(loss, working_model), psi, design_matrix, Q)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  H0 <- torch_tensor(-1 / (1 - nuisance$pi))
  H1 <- torch_tensor(1 / nuisance$pi)
  H  <- torch_tensor(ifelse(data[[A]] == 1, 1 / nuisance$pi, -1 / (1 - nuisance$pi)))

  ##### One-step estimator
  onestep_est <- onestep(
    Lm(loss, working_model), psi, beta, design_matrix, Q, Delta(H, Yt, torch::torch_tensor(nuisance$mu))
  )

  draws <- 1e3
  onestep_joint <- mvtnorm::rmvnorm(draws, mean = as.numeric(onestep_est$est), sigma = var(as.matrix(onestep_est$eif)) / n)

  ##### TMLE
  if(tmle == TRUE) {
    # Calculate clever covariates for fluctuation model
    calculate_clever <- \(Lm, H, H0, H1, psi, Q, beta, design_matrix) {
      p <- ncol(design_matrix)
      n <- nrow(design_matrix)

      clever <- torch::torch_tensor(matrix(0, n, p))
      clever0 <- torch::torch_tensor(matrix(0, n, p))
      clever1 <- torch::torch_tensor(matrix(0, n, p))
      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      for(i in 1:n) {
        d <- Minv$matmul(grad_dL(Lm, psi[i], beta, design_matrix[i, ]))
        clever[i, ]  <- torch::torch_reshape(H[i] * d, p)
        clever0[i, ] <- torch::torch_reshape(H0[i] * d, p)
        clever1[i, ] <- torch::torch_reshape(H1[i] * d, p)
      }
      list(
        clever = clever,
        clever0 = clever0,
        clever1 = clever1
      )
    }

    # Clever covariate for marginal distribution of covariates
    calculate_K <- \(Lm, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]

      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      K <- torch::torch_tensor(matrix(0, n, p))
      for(i in 1:n) {
        K[i, ] <- Minv$matmul(dL(Lm, psi[i], beta, design_matrix[i, ])[[1]])$reshape(p)
      }
      K
    }

    dQ_fluctuation_depsilon <- \(epsilon, K, Q) {
      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      Q_fluctuation(epsilon, K, Q)$reshape(c(n, 1))$mul(K - Q$reshape(c(n, 1))$mul(K)$mul(exp(K * epsilon)) / Qn)
    }

    Q_fluctuation <- \(epsilon, K, Q) {
      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      Q  <- exp(K$matmul(epsilon) * Q) / Qn
      Q
    }

    # TMLE fluctuation model
    if(tmle_linear == TRUE) {
      tmle_loss <- nn_mse_loss(reduction = "sum")
    }
    else {
      tmle_loss <- nn_bce_with_logits_loss(reduction = "sum")
    }

    tmle_fluctuation_model <- \(epsilon, mu, mu0, mu1, clever, clever0, clever1, K, Q, Y, design_matrix, Lm = NULL, condvar = NA, bayes = FALSE) {
      # Fluctuate covariate distribution
      Q <- Q_fluctuation(epsilon, K, Q)

      # Fluctuate outcome regression
      if(tmle_linear == TRUE) {
        mu  <- mu  + clever$matmul(epsilon)
        mu0 <- mu0 + clever0$matmul(epsilon)
        mu1 <- mu1 + clever1$matmul(epsilon)
        target <- tmle_loss(mu, Y)
      }
      else {
        mu_logit <- mu$logit() + clever$matmul(epsilon)
        mu0      <- torch::torch_sigmoid(mu0$logit() + clever0$matmul(epsilon))
        mu1      <- torch::torch_sigmoid(mu1$logit() + clever1$matmul(epsilon))
        target   <- tmle_loss(mu_logit, Y)
      }

      # Combined loss function
      target <- target - as.numeric(log(Q)$sum())

      if(bayes == FALSE) {
        return(target)
      }
      else {
        psi <- mu1 - mu0
        beta <- B(Lm(loss, working_model), psi$detach()$clone(), design_matrix, Q$detach()$clone())

        if(tmle_linear == TRUE) {
          target <- as.numeric(-((mu - Y)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()) + as.numeric(log(Q)$sum())
        }
        else {
          target <- -as.numeric(target)
        }

        # Prior
        target <- target + bayes_prior(as.numeric(beta))
        # Jacobian adjustment
        jacobian <- torch::torch_transpose(dB_dpsi(Lm(loss, working_model), psi, Q, design_matrix, beta), 1, 2)
        if(tmle_linear == TRUE) {
          jacobian <- jacobian$matmul(clever1 - clever0)
        }
        else {
          jacobian <- jacobian$matmul(clever1 * mu1$reshape(c(n, 1)) * (1 - mu1$reshape(c(n, 1))) - clever0 * mu0$reshape(c(n, 1)) * (1 - mu0$reshape(c(n, 1))))
        }
        jacobian <- jacobian + torch::torch_transpose(dB_dQ(Lm(loss, working_model), psi, Q, design_matrix, beta), 1, 2)$matmul(dQ_fluctuation_depsilon(epsilon, K, Q))
        target <- target + log(abs(jacobian$det()))

        ret <- list(log.density = as.numeric(target), beta = as.numeric(beta))

        return(ret)
      }
    }

    Q_star <- Q
    mu_star  <- torch::torch_tensor(nuisance$mu)
    mu0_star <- torch::torch_tensor(nuisance$mu0, requires_grad = TRUE)
    mu1_star <- torch::torch_tensor(nuisance$mu1, requires_grad = TRUE)
    beta_star <- beta
    epsilon_star <- rep(0, p)

    converged <- TRUE

    for(tmle_iter in 1:tmle_maxiter) {
      psi_star <- (mu1_star - mu0_star)$detach()$clone()
      psi_star$requires_grad_(TRUE)

      K <- calculate_K(Lm(loss, working_model), psi_star, Q_star, beta_star, design_matrix)
      clever <- calculate_clever(Lm(loss, working_model), H, H0, H1, psi_star, Q_star, beta_star, design_matrix)
      epsilon_star <- tmle_mle(p, tmle_fluctuation_model, mu_star, mu0_star, mu1_star, clever$clever, clever$clever0, clever$clever1, K, Q_star, Yt, design_matrix)

      if(any(is.nan(as.numeric(epsilon_star)))) {
        warning("TMLE failed to converge.")
        converged <- FALSE
        break
      }

      m <- max(as.numeric(epsilon_star))
      cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      if(tmle_linear == TRUE) {
        mu_star  <- mu_star  + clever$clever$matmul(epsilon_star)
        mu0_star <- mu0_star + clever$clever0$matmul(epsilon_star)
        mu1_star <- mu1_star + clever$clever1$matmul(epsilon_star)
      }
      else {
        mu_star  <- torch::torch_sigmoid(mu_star$logit()  + clever$clever$matmul(epsilon_star))
        mu0_star <- torch::torch_sigmoid(mu0_star$logit() + clever$clever0$matmul(epsilon_star))
        mu1_star <- torch::torch_sigmoid(mu1_star$logit() + clever$clever1$matmul(epsilon_star))
      }

      Q_star  <- Q_fluctuation(epsilon_star, K, Q_star)

      beta_star <- B(Lm(loss, working_model), psi_star, design_matrix, Q_star)$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if(abs(m) < 1e-3) {
        break
      }
    }

    tmle_est <- tmle_se <- tmle_lower <- tmle_upper <- rep(NA, p)
    tmle_eif <- matrix(NA, ncol = p, nrow = n)
    if(converged == TRUE) {
      tmle_est   <- B(Lm(loss, working_model), psi_star, design_matrix, Q_star)
      tmle_eif   <- eif(Lm(loss, working_model), psi_star, tmle_est, design_matrix, Q_star, Delta(H, Yt, mu_star))
      tmle_se    <- apply(tmle_eif, 2, sd) / sqrt(n)
      tmle_lower <- tmle_est + qnorm(0.025) * tmle_se
      tmle_upper <- tmle_est + qnorm(0.975) * tmle_se
      tmle_joint <- mvtnorm::rmvnorm(draws, mean = as.numeric(as.matrix(tmle_est)), sigma = var(as.matrix(tmle_eif)) / n)
    }

    tmle_beta_samples <- NULL
    tmle_acc_rate <- NULL

    if(bayes == TRUE) {
      tmle_beta_samples <- array(dim = c(4, bayes_draws, p))
      tmle_acc_rate <- 0

      for(chain in 1:4) {
        cat("Chain: ", chain, "\n\n")

        log_dens <- \(epsilon) tmle_fluctuation_model(
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

        mcmc <- adaptMCMC::MCMC(
          log_dens,
          n = bayes_draws,
          init = as.numeric(epsilon_star),
          adapt = TRUE,
          acc.rate = 0.3,
          scale = rep(1e-3, p)
        )

        tmle_beta_samples[chain, ,] <- matrix(unlist(mcmc$extra.values), ncol = p, nrow = bayes_draws, byrow = TRUE)
        tmle_acc_rate <- tmle_acc_rate + 1/4 * mcmc$acceptance.rate
      }

    }
  }

  res <- list(
    estimand = "treatment_effect_modification",
    p = p,
    f = f,
    terms = terms,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    nuisance = nuisance,
    plugin = list(
      est = as.numeric(plugin)
    ),
    onestep = list(
      est   = as.numeric(onestep_est$est),
      se    = as.numeric(onestep_est$se),
      lower = as.numeric(onestep_est$lower),
      upper = as.numeric(onestep_est$upper),
      eif   = onestep_est$eif,
      psi   = as.numeric(psi),
      joint_draws = onestep_joint
    )
  )

  if(tmle == TRUE) {
    res$tmle <- list(
      est   = as.numeric(tmle_est),
      se    = as.numeric(tmle_se),
      lower = as.numeric(tmle_lower),
      upper = as.numeric(tmle_upper),
      eif   = tmle_eif,
      psi   = as.numeric(psi_star),
      joint_draws = tmle_joint,
      samples = tmle_beta_samples,
      acc_rate = tmle_acc_rate
    )
  }

  class(res) <- "targeted_msm"

  res
}
