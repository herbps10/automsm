estimate_categorical_dose_response_nuisance <- function(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type) {
  n <- nrow(data)
  As <- sort(unique(data[[A]]))
  k <- length(As)
  pi_a_hat <- mu_a_hat <- matrix(0, n, k)
  pi_hat <- mu_hat <- numeric(n)

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)
  outcome_family <- stats::gaussian()
  if(outcome_type == "binomial") outcome_family <- stats::binomial()

  if(outer_folds > 1) {
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
    }
  }

  for(a_index in seq_along(As)) {
    ind <- data[[A]] == As[a_index]
    mu_hat[ind] <- mu_a_hat[ind, a_index]
    pi_hat[ind] <- pi_a_hat[ind, a_index]
  }

  list(
    pi_a = pi_a_hat,
    mu_a = mu_a_hat,
    mu   = mu_hat,
    pi   = pi_hat
  )
}

categorical_dose_response <- function(data, X, A, Y, formula, loss = loss_squared_error, working_model = working_model_linear, learners_trt = "glm", learners_outcome = "glm", outer_folds = 5, inner_folds = 5, tmle = TRUE, tmle_maxiter = 25, outcome_type = "binomial") {
  n <- nrow(data)
  #data <- data[, c(X, A, Y)]
  nuisance <- estimate_categorical_dose_response_nuisance(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, outcome_type)

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
  for(k in 1:K) {
    H[, k] <- (data[[A]] == As[k]) / nuisance$pi
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
      onestep_model[index, k] <- as.numeric(working_model(params[index, ], design_matrix[k, 1, ]))
    }
  }

  ##### TMLE
  if(tmle == TRUE) {
    # Calculate clever covariates for fluctuation model
    calculate_clever <- \(Lm, H, psi, Q, beta, design_matrix) {
      p <- rev(dim(design_matrix))[1]
      n <- rev(dim(design_matrix))[2]

      clever <- torch::torch_tensor(matrix(0, n, p))
      Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)
      for(i in 1:n) {
        d <- Minv$matmul(grad_dL(Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE]))
        clever[i, ]  <- torch::torch_reshape(d$matmul(H[i,]), p)
      }
      list(
        clever = clever
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

    # TMLE fluctuation model
    mse_loss <- nn_mse_loss(reduction = "sum")
    tmle_fluctuation_model <- \(epsilon, mu, clever, K, Q, Y) {
      mu  <- mu  + clever$matmul(epsilon)

      Qn <- exp(K$matmul(epsilon) * Q)$sum()
      Q  <- exp(K$matmul(epsilon) * Q) / Qn

      target <- mse_loss(mu, Y)
      target <- target + log(Q)$sum()
      target
    }

    Q_star <- Q
    mu_star  <- torch::torch_tensor(nuisance$mu)
    mu_a_star  <- torch::torch_tensor(nuisance$mu_a)
    beta_star <- beta
    epsilon_star <- rep(0, p)

    for(tmle_iter in 1:tmle_maxiter) {
      psi_star <- mu_a_star$detach()$clone()
      psi_star$requires_grad_(TRUE)

      clever_K <- calculate_K(combined_loss, psi_star, Q_star, beta_star, design_matrix)
      clever <- calculate_clever(combined_loss, H, psi_star, Q_star, beta_star, design_matrix)
      epsilon_star <- tmle_mle(p, tmle_fluctuation_model, mu_star, clever$clever, clever_K, Q_star, Yt)

      m <- max(as.numeric(epsilon_star))
      cat(glue::glue("TMLE iteration: {tmle_iter}, max(epsilon): {m}\n\n"))

      mu_star  <- mu_star + clever$clever$matmul(epsilon_star)
      for(k in 1:K) {
        mu_a_star[, k] <- mu_a_star[, k] + clever$clever$matmul(epsilon_star)
      }

      Qn <- exp(clever_K$matmul(epsilon_star) * Q_star)$sum()
      Q_star  <- exp(clever_K$matmul(epsilon_star) * Q_star) / Qn

      beta_star <- B(combined_loss, psi_star, design_matrix, Q_star)$detach()$clone()
      beta_star$requires_grad_(TRUE)

      if(abs(m) < 1e-3) {
        break
      }
    }

    tmle_est   <- B(combined_loss, psi_star, design_matrix, Q_star)
    tmle_eif   <- eif(combined_loss, psi_star, tmle_est, design_matrix, Q_star, Delta(H, Yt, mu_star))
    tmle_se    <- apply(tmle_eif, 2, sd) / sqrt(n)
    tmle_lower <- tmle_est + qnorm(0.025) * tmle_se
    tmle_upper <- tmle_est + qnorm(0.975) * tmle_se
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
      psi   = as.numeric(psi_star),
      model = as.numeric(working_model(tmle_est, design_matrix))
    )
  }

  class(res) <- "targeted_msm"

  res
}
