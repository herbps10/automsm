msm_spec_dose_response <- function(tmle_linear = TRUE) {
  tmle_loss <- if(tmle_linear) {
    torch::nn_mse_loss(reduction = "sum")
  }
  else {
    torch::nn_bce_with_logits_loss(reduction = "sum")
  }

  new_msm_spec(
    estimand = "dose_response",
    tol = 1e-2,
    nan_guard = TRUE,
    supports_bayes = TRUE,
    tmle_loss = tmle_loss,

    sweep_criterion = function(eps_list) max(abs(as.numeric(eps_list[[1]]))),

    init_state = function(problem) {
      list(
        mu = list(
          obs = as_float_tensor(problem$nuisance$mu),
          arms = as_float_tensor(problem$nuisance$mu_a)
        ),
        Q = problem$Q0,
        beta = NULL
      )
    },

    steps = function(problem) list(list(id = 1L, fluctuate_Q = TRUE)),

    psi_from_state = function(problem, state, detach = TRUE) {
      psi <- state$mu$arms
      if(detach) {
        psi <- psi$detach()$clone()
        psi$requires_grad_(TRUE)
      }
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      K <- problem$K
      p <- problem$p
      H <- problem$aux$H
      HA <- problem$aux$HA
      cl <- torch::torch_zeros(c(n, p))
      clA <- torch::torch_zeros(c(n, K, p))
      for(i in 1:n) {
        d <- Minv$matmul(grad_dL(
          problem$Lm_fn, psi[i, drop = FALSE], state$beta,
          problem$design_matrix[i, , drop = FALSE]
        ))

        cl[i, ] <- torch::torch_reshape(d$matmul(H[i, ]), p)
        for(j in 1:K) {
          x <- torch::torch_zeros(K)
          x[j] <- HA[i, j]
          clA[i, j, ] <- torch::torch_reshape(d$matmul(x), p)
        }
      }
      list(obs = cl, arms = clA)
    },

    fluctuation_objective = function(epsilon, problem, step, state, clever, K_Q) {
      pred <- if(tmle_linear) {
        state$mu$obs + clever$obs$matmul(epsilon)
      }
      else {
        state$mu$obs$logit() + clever$obs$matmul(epsilon)
      }
      Qf <- Q_fluctuation(epsilon, K_Q, state$Q)
      tmle_loss(pred, problem$Yt) - Qf$log()$sum()
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q) {
      K <- problem$K
      obs <- if(tmle_linear) {
        state$mu$obs + clever$obs$matmul(epsilon)
      }
      else {
        torch::torch_sigmoid(state$mu$obs$logit() + clever$obs$matmul(epsilon))
      }

      arms <- state$mu$arms$detach()$clone()
      for(k in 1:K) {
        arms[, k] <- if(tmle_linear) {
          state$mu$arms[, k] + clever$arms[, k, ]$matmul(epsilon)
        }
        else {
          torch::torch_sigmoid(
            state$mu$arms[, k]$logit() + clever$arms[, k, ]$matmul(epsilon)
          )
        }
      }
      state$mu <- list(obs = obs$detach(), arms = arms$detach())
      if(isTRUE(step$fluctuate_Q)) {
        state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
      }
      state
    },

    delta = function(problem, state) {
      problem$aux$H$mul(
        (problem$Yt - state$mu$obs)$reshape(c(problem$n, 1))
      )
    },

    dpsi_depsilon = function(problem, state, clever, epsilon) {
      n <- problem$n
      lapply(seq_len(problem$K), function(k) {
        if(tmle_linear) {
          clever$arms[, k, ]
        }
        else {
          m <- state$mu$arms[, k]$reshape(c(n, 1))
          clever$arms[, k, ] * m * (1 - m)
        }
      })
    },

    bayes_loglik = function(epsilon, problem, state, clever, Q_eps, condvar) {
      if(tmle_linear) {
        pred <- state$mu$obs + clever$obs$matmul(epsilon)
        stopifnot(!is.null(condvar))
        as.numeric(-((pred - problem$Yt)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum()) + as.numeric(Q_eps$log()$sum())
      }
      else {
        lg <- state$mu$obs$logit() + clever$obs$matmul(epsilon)
        -as.numeric(tmle_loss(lg, problem$Yt)) + as.numeric(Q_eps$log()$sum())
      }
    }
  )
}
