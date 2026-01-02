estimate_treatment_effect_modification_nuisance <- function(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type) {
  n <- nrow(data)
  data0 <- data1 <- data
  data0[[A]] <- 0
  data1[[A]] <- 1
  pi_hat <- mu0_hat <- mu1_hat <- numeric(n)

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)
  outcome_family <- stats::gaussian()
  if(outcome_type == "binomial") outcome_family <- stats::binomial()

  if(outer_folds > 1) {
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
    }
  }
  else {
    pi_model <- SuperLearner::SuperLearner(
      Y = data[[A]],
      X = data[, X, drop = FALSE],
      SL.library = learners_trt,
      cvControl = cv_control,
      family = "binomial",
      env = environment(SuperLearner::SuperLearner)
    )

    mu_model <- SuperLearner::SuperLearner(
      Y = data[[Y]],
      X = data[, c(X, A), drop = FALSE],
      SL.library = learners_outcome,
      family = outcome_family,
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    pi_hat  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data, onlySL = TRUE)$pred
    mu0_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data0, onlySL = TRUE)$pred
    mu1_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data1, onlySL = TRUE)$pred
  }
  mu_hat <- ifelse(data[[A]] == 1, mu1_hat, mu0_hat)

  list(
    pi  = pi_hat,
    mu0 = mu0_hat,
    mu1 = mu1_hat,
    mu  = mu_hat
  )
}

treatment_effect_modification <- function(data, X, A, Y, formula, loss = loss_squared_error, working_model = working_model_linear, learners_trt = "glm", learners_outcome = "glm", outer_folds = 5, inner_folds = 5, tmle = TRUE, tmle_maxiter = 25, outcome_type = "binomial") {
  n <- nrow(data)
  #data <- data[, c(X, A, Y)]
  nuisance <- estimate_treatment_effect_modification_nuisance(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type)

  Xt <- torch::torch_tensor(as.matrix(data[, X, drop = FALSE]))
  Yt <- torch::torch_tensor(data[[Y]])

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

    # TMLE fluctuation model
    mse_loss <- nn_mse_loss(reduction = "sum")
    tmle_fluctuation_model <- \(epsilon, mu, mu0, mu1, clever, clever0, clever1, K, Q, Y) {
      mu  <- mu  + clever$matmul(epsilon)

      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      Q  <- exp(K$matmul(epsilon) * Q) / Qn

      #target <- (mu - Y)$pow(2)$sum()
      target <- mse_loss(mu, Y)
      target <- target + log(Q)$sum()
      target
    }

    Q_star <- Q
    mu_star  <- torch::torch_tensor(nuisance$mu)
    mu0_star <- torch::torch_tensor(nuisance$mu0, requires_grad = TRUE)
    mu1_star <- torch::torch_tensor(nuisance$mu1, requires_grad = TRUE)
    beta_star <- beta
    epsilon_star <- rep(0, p)

    for(tmle_iter in 1:tmle_maxiter) {
      psi_star <- (mu1_star - mu0_star)$detach()$clone()
      psi_star$requires_grad_(TRUE)

      K <- calculate_K(Lm(loss, working_model), psi_star, Q_star, beta_star, design_matrix)
      clever <- calculate_clever(Lm(loss, working_model), H, H0, H1, psi_star, Q_star, beta_star, design_matrix)
      epsilon_star <- tmle_mle(p, tmle_fluctuation_model, mu_star, mu0_star, mu1_star, clever$clever, clever$clever0, clever$clever1, K, Q_star, Yt)

      m <- max(as.numeric(epsilon_star))
      cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      mu_star  <- mu_star  + clever$clever$matmul(epsilon_star)
      mu0_star <- mu0_star + clever$clever0$matmul(epsilon_star)
      mu1_star <- mu1_star + clever$clever1$matmul(epsilon_star)

      Qn <- exp(K$matmul(epsilon_star) * Q_star)$sum()
      Q_star  <- exp(K$matmul(epsilon_star) * Q_star) / Qn

      beta_star <- B(Lm(loss, working_model), psi_star, design_matrix, Q_star)$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if(abs(m) < 1e-3) {
        break
      }
    }

    tmle_est   <- B(Lm(loss, working_model), psi_star, design_matrix, Q_star)
    tmle_eif   <- eif(Lm(loss, working_model), psi_star, tmle_est, design_matrix, Q_star, Delta(H, Yt, mu_star))
    tmle_se    <- apply(tmle_eif, 2, sd) / sqrt(n)
    tmle_lower <- tmle_est + qnorm(0.025) * tmle_se
    tmle_upper <- tmle_est + qnorm(0.975) * tmle_se
  }

  res <- list(
    estimand = "treatment_effect_modification",
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
      model = as.numeric(working_model(onestep_est$est, design_matrix))
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
      model = as.numeric(working_model(tmle_est, design_matrix))
    )
  }

  class(res) <- "targeted_msm"

  res
}
