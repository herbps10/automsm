estimate_categorical_dose_response_nuisance <- function(
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
  if(outer_folds == 1) cv[[1]]$training_set <- cv[[1]]$validation_set
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)
  outcome_family <- stats::gaussian()
  if(outcome_type == "binomial") outcome_family <- stats::binomial()

  for(fold in seq_along(cv)) {
    training   <- cv[[fold]]$training_set
    validation <- cv[[fold]]$validation_set

    for(a_index in seq_along(As)) {
      pi_model <- SuperLearner::SuperLearner(
        Y = as.numeric(data[[A]][training] == As[a_index]),
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      pi_a_hat[validation, a_index] <- SuperLearner::predict.SuperLearner(pi_model, newdata = data[validation, c(X), drop = FALSE], onlySL = TRUE)$pred
    }

    mu_model <- SuperLearner::SuperLearner(
      Y = data[[Y]][training],
      X = data[training, c(X, A), drop = FALSE],
      SL.library = learners_outcome,
      family = outcome_family,
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    for(a_index in seq_along(As)) {
      newdata <- data[validation, c(X, A)]
      newdata[[A]] <- As[a_index]
      mu_a_hat[validation, a_index] <- SuperLearner::predict.SuperLearner(mu_model, newdata = newdata, onlySL = TRUE)$pred
    }

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
      condvar_hat[validation] <- ifelse(condvar_hat[validation] < epsilon, epsilon, condvar_hat[validation])
    }
  }

  for(a_index in seq_along(As)) {
    ind <- data[[A]] == As[a_index]
    mu_hat[ind] <- mu_a_hat[ind, a_index]
    pi_hat[ind] <- pi_a_hat[ind, a_index]
  }

  list(
    pi_a    = pi_a_hat,
    mu_a    = mu_a_hat,
    mu      = mu_hat,
    pi      = pi_hat,
    condvar = condvar_hat
  )
}

#' @export
categorical_dose_response <- function(
    data,
    X,
    A,
    Y,
    formula,
    loss = loss_squared_error,
    working_model = working_model_linear,
    learners_trt = "glm",
    learners_outcome = "glm",
    outer_folds = 5,
    inner_folds = 5,
    tmle = TRUE,
    tmle_maxiter = 25,
    outcome_type = "binomial",
    bayes = FALSE,
    bayes_draws = 1e3,
    bayes_prior = \(beta) sum(dnorm(as.numeric(beta), mean = 0, sd = 1, log = TRUE)),
    epsilon = 1e-5,
    nuisance = NULL
) {
  n <- nrow(data)
  #data <- data[, c(X, A, Y)]
  if(is.null(nuisance)) {
    nuisance <- estimate_categorical_dose_response_nuisance(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type, estimate_conditional_variance = bayes, epsilon = epsilon)
  }

  Xt <- torch::torch_tensor(as.matrix(data[, X, drop = FALSE]))
  Yt <- torch::torch_tensor(data[[Y]])

  As <- sort(unique(data[[A]]))
  K <- length(As)

  mat <- model.matrix(formula, data = data)
  terms <- colnames(mat)
  p <- ncol(mat)

  design_matrix <- torch::torch_zeros(K * n * ncol(mat))$reshape(c(K, n, ncol(mat)))

  for(k in 1:K) {
    datak <- data
    datak[[A]] <- As[k]
    mat <- model.matrix(formula, data = datak)
    terms <- colnames(mat)
    design_matrix[k,,] <- mat
  }

  combined_loss <- \(t, beta, X) {
    n <- rev(dim(X))[2]
    K <- dim(X)[1]
    sum <- torch::torch_tensor(rep(0, n), requires_grad = TRUE)
    for(k in 1:K) {
      sum <- sum$add(loss(t[, k], working_model(beta, X[k,,])))
    }
    sum
  }

  Delta <- \(H, Y, mu) {
    n <- nrow(H)
    H$mul((Y - mu)$reshape(c(n, 1)))
  }

  ##### Plugin estimator
  Q <- torch::torch_tensor(rep(1 / n, n))
  psi <- torch::torch_tensor(nuisance$mu_a, requires_grad = TRUE)
  plugin <- B(combined_loss, psi, design_matrix, Q)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  H <- torch_zeros(c(n, K))
  HA <- torch_zeros(c(n, K))
  for(k in 1:K) {
    H[, k] <- (data[[A]] == As[k]) / nuisance$pi
    HA[, k] <- 1 / nuisance$pi_a[, k]
  }

  ##### One-step estimator
  onestep_est <- onestep(
    combined_loss, psi, beta, design_matrix, Q, Delta(H, Yt, torch::torch_tensor(nuisance$mu))
  )


  draws <- 1e3
  onestep_model <- matrix(nrow = draws, ncol = K)
  params <- mvtnorm::rmvnorm(draws, mean = as.numeric(onestep_est$est), sigma = var(as.matrix(onestep_est$eif)))
  for(index in 1:draws) {
    for(k in 1:K) {
      onestep_model[index, k] <- as.numeric(working_model(torch::torch_tensor(params[index, ]), design_matrix[k, 1, ]$reshape(c(1, p))))
    }
  }

  ##### TMLE
  if(tmle == TRUE) {
    # Calculate clever covariates for fluctuation model
    calculate_clever <- \(Lm, H, HA, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]

      clever <- torch::torch_tensor(matrix(0, n, p))
      cleverA <- torch::torch_tensor(array(0, dim = c(K, n, p)))
      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      for(i in 1:n) {
        d <- Minv$matmul(grad_dL(Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE]))
        clever[i, ]  <- torch::torch_reshape(d$matmul(H[i,]), p)
        for(k in 1:K) {
          x <- torch::torch_zeros(K)
          x[k] <- HA[i, k]
          cleverA[k, i, ] <- torch::torch_reshape(d$matmul(x), p)
        }
      }
      list(
        clever = clever,
        cleverA = cleverA
      )
    }

    # Clever covariate for marginal distribution of covariates
    calculate_K <- \(Lm, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]

      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      K <- torch::torch_tensor(matrix(0, n, p))
      for(i in 1:n) {
        K[i, ] <- Minv$matmul(dL(Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE])[[1]])$reshape(p)
      }
      K
    }

    dQ_fluctuation_depsilon <- \(epsilon, K, Q) {
      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      Q_fluctuation(epsilon, K, Q)$reshape(c(n, 1))$mul(K - Q$reshape(c(n, 1))$mul(K)$mul(exp(K * epsilon)) / Qn)
    }

    Q_fluctuation <- \(epsilon, K, Q) {
      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      exp(K$matmul(epsilon) * Q) / Qn
    }

    # TMLE fluctuation model
    tmle_linear <- TRUE
    if(tmle_linear == TRUE) {
      tmle_loss <- nn_mse_loss(reduction = "sum")
    }
    else {
      tmle_loss <- nn_bce_with_logits_loss(reduction = "sum")
    }

    tmle_fluctuation_model <- \(epsilon, mu, mu_a, clever, clever_K, Q, Y, design_matrix, condvar = NULL, bayes = FALSE) {
      mu  <- mu + clever$clever$matmul(epsilon)
      Q <- Q_fluctuation(epsilon, clever_K, Q)

      if(bayes == FALSE) {
        target <- tmle_loss(mu, Y)
        target <- target - log(Q)$sum()
        target
      }
      else {
        mu_a <- mu_a + clever$cleverA$matmul(epsilon)$transpose(1, 2)

        beta <- B(combined_loss, mu_a, design_matrix, Q)

        target <- as.numeric(-((mu - Y)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()) + as.numeric(log(Q)$sum())

        # Prior
        #target <- target + bayes_prior(as.numeric(beta))
        #jacobian <- torch::torch_zeros(c(p, p))
        #J <- dB_dpsi(combined_loss, mu_a, Q, design_matrix, beta)
        #for(k in 1:K) {
        #  jacobian <- jacobian + J[k, ,]$transpose(1, 2)$matmul(clever$cleverA[k,])
        #}
        #jacobian <- jacobian + torch::torch_transpose(dB_dQ(combined_loss, psi, Q, design_matrix, beta), 1, 2)$matmul(dQ_fluctuation_depsilon(epsilon, clever_K, Q))
        #target <- target + log(abs(jacobian$det()))

        ret <- list(log.density = as.numeric(target), beta = as.numeric(beta))

        return(ret)
      }
    }

    Q_star <- Q
    mu_star  <- torch::torch_tensor(nuisance$mu)
    mu_a_star  <- torch::torch_tensor(nuisance$mu_a)
    beta_star <- beta
    epsilon_star <- rep(0, p)

    for(tmle_iter in 1:tmle_maxiter) {
      mu_a_star <- mu_a_star$detach()$clone()
      mu_a_star$requires_grad_(TRUE)

      clever_K <- calculate_K(combined_loss, mu_a_star, Q_star, beta_star, design_matrix)
      clever <- calculate_clever(combined_loss, H, HA, mu_a_star, Q_star, beta_star, design_matrix)
      epsilon_star <- tmle_mle(p, tmle_fluctuation_model, mu_star, mu_a_star, clever, clever_K, Q_star, Yt)

      m <- max(as.numeric(epsilon_star))
      cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      mu_star <- mu_star + clever$clever$matmul(epsilon_star)
      mu_a_star <- mu_a_star$detach()$clone()
      for(k in 1:K) {
        mu_a_star[, k] <- mu_a_star[, k] + clever$cleverA[k, ,]$matmul(epsilon_star)
      }

      Q_star <- Q_fluctuation(epsilon_star, clever_K, Q_star)

      beta_star <- B(combined_loss, mu_a_star$detach(), design_matrix, Q_star$detach())$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if(abs(m) < 1e-2) {
        break
      }
    }

    mu_a_star <- mu_a_star$detach()$clone()
    mu_a_star$requires_grad_(TRUE)

    tmle_est   <- B(combined_loss, mu_a_star$detach(), design_matrix, Q_star$detach())
    tmle_eif   <- eif(combined_loss, mu_a_star, tmle_est, design_matrix, Q_star$detach(), Delta(H, Yt, mu_star$detach()))
    tmle_se    <- apply(tmle_eif, 2, sd) / sqrt(n)
    tmle_lower <- tmle_est + qnorm(0.025) * tmle_se
    tmle_upper <- tmle_est + qnorm(0.975) * tmle_se

    tmle_beta_samples <- NULL
    tmle_acc_rate <- NULL

    if(bayes == TRUE) {
      tmle_beta_samples <- array(dim = c(4, bayes_draws, p))
      tmle_acc_rate <- 0

      for(chain in 1:4) {
        log_dens <- \(epsilon) tmle_fluctuation_model(
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
    estimand = "categorical_dose_response",
    p = p,
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
      model = onestep_model
    )
  )

  if(tmle == TRUE) {
    res$tmle <- list(
      est   = as.numeric(tmle_est),
      se    = as.numeric(tmle_se),
      lower = as.numeric(tmle_lower),
      upper = as.numeric(tmle_upper),
      eif   = tmle_eif,
      samples = tmle_beta_samples,
      acc_rate = tmle_acc_rate
    )
  }

  class(res) <- "targeted_msm"

  res
}
