msm_spec_cate <- function(tmle_linear = TRUE) {
  tmle_loss <- if(tmle_linear) {
    torch::nn_mse_loss(reduction = "sum")
  }
  else {
    torch::nn_bce_with_logits_loss(reduction = "sum")
  }

  fluctuate <- function(mu_col, cc, epsilon) {
    if(tmle_linear) {
      mu_col + cc$matmul(epsilon)
    }
    else {
      torch::torch_sigmoid(mu_col$logit() + cc$matmul(epsilon))
    }
  }

  new_msm_spec(
    estimand = "cate",
    tol = 1e-3,
    nan_guard = TRUE,
    supports_bayes = TRUE,
    tmle_loss = tmle_loss,

    sweep_criterion = function(eps_list) max(abs(as.numeric(eps_list[[1]]))),

    init_state = function(problem) {
      n <- problem$n
      list(
        mu = list(
          obs = as_float_tensor(problem$nuisance$mu)$reshape(c(n, 1)),
          a0 = as_float_tensor(problem$nuisance$mu0)$reshape(c(n, 1)),
          a1 = as_float_tensor(problem$nuisance$mu1)$reshape(c(n, 1))
        ),
        Q = problem$Q0,
        beta = NULL
      )
    },

    steps = function(problem) list(list(id = 1L, fluctuate_Q = TRUE)),

    psi_from_state = function(problem, state, detach = TRUE) {
      psi <- state$mu$a1 - state$mu$a0
      if(detach) {
        psi <- psi$detach()$clone()
        psi$requires_grad_(TRUE)
      }
      psi
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      p <- problem$p
      H <- problem$aux$H
      H0 <- problem$aux$H0
      H1 <- problem$aux$H1
      cl <- torch::torch_zeros(c(n, p))
      cl0 <- torch::torch_zeros(c(n, p))
      cl1 <- torch::torch_zeros(c(n, p))

      for(i in 1:n) {
        d <- Minv$matmul(grad_dL(
          problem$Lm_fn, psi[i, drop = FALSE], state$beta,
          problem$design_matrix[i, , drop = FALSE]
        ))

        cl[i, ] <- torch::torch_reshape(H[i] * d, p)
        cl0[i, ] <- torch::torch_reshape(H0[i] * d, p)
        cl1[i, ] <- torch::torch_reshape(H1[i] * d, p)
      }
      list(obs = cl, a0 = cl0, a1 = cl1)
    },

    fluctuation_objective = function(epsilon, problem, step, state, clever, K_Q) {
      if(tmle_linear) {
        pred <- state$mu$obs[, 1] + clever$obs$matmul(epsilon)
      }
      else {
        pred <- state$mu$obs[, 1]$logit() + clever$obs$matmul(epsilon)
      }
      target <- tmle_loss(pred, problem$Yt)
      target - state$Q$log()$sum()
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q) {
      n <- problem$n
      state$mu <- list(
        obs = fluctuate(state$mu$obs[, 1], clever$obs, epsilon)$reshape(c(n, 1)),
        a0 = fluctuate(state$mu$a0[, 1], clever$a0, epsilon)$reshape(c(n, 1)),
        a1 = fluctuate(state$mu$a1[, 1], clever$a1, epsilon)$reshape(c(n, 1))
      )

      if(isTRUE(step$fluctuate_Q)) {
        state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
      }

      state
    },

    delta = function(problem, state) {
      (problem$aux$H * (problem$Yt - state$mu$obs[, 1]))$reshape(c(problem$n, 1))
    },

    dpsi_depsilon <- function(problem, state, clever, epsilon) {
      n <- problem$n
      if(tmle_linear) {
        list(clever$a1 - clever$a0)
      }
      else {
        m1 <- state$mu$a1$reshape(c(n, 1))
        m0 <- state$mu$a0$reshape(c(n, 1))
        list(clever$a1 * m1 * (1 - m1) - clever$a0 * m0 * (1 - m0))
      }
    },

    bayes_loglik = function(epsilon, problem, state, clever, Q_eps, condvar) {
      pred <- state$mu$obs[, 1] + clever$obs$matmul(epsilon)

      if(tmle_linear) {
        stopifnot(!is.null(condvar))
        as.numeric(-((pred - problem$Yt)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()) + as.numeric(Q_eps$log()$sum())
      }
      else {
        lg <- state$mu$obs[, ]$logit() + clever$obs$matmul(epsilon)
        -as.numeric(tmle_loss(lg, problem$Yt)) + as.numeric(Q_eps$log()$sum())
      }
    }
  )
}
