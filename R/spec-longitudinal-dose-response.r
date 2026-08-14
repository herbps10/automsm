msm_spec_longitudinal_dose_response <- function(tmle_linear = TRUE) {
  tmle_loss <- if(tmle_linear) {
    torch::nn_mse_loss(reduction = "sum")
  }
  else {
    torch::nn_bce_with_logits_loss(reduction = "sum")
  }

  new_msm_spec(
    estimand = "longitudinal_dose_response",
    tol = 1e-2,
    nan_guard = TRUE,
    supports_bayes = FALSE,
    tmle_loss = tmle_loss,

    sweep_criterion = function(eps_list) {
      max(vapply(eps_list, function(e) max(abs(as.numeric(e))), numeric(1)))
    },

    init_state = function(problem) {
      list(
        mu = list(
          nodes = lapply(seq_len(problem$tau + 1L), function(t) {
            as_float_tensor(t(problem$nuisance$mu[, , t]))
          })
        ),
        Q = problem$Q0,
        beta = NULL,
        initial = TRUE
      )
    },

    # Backward sweep over ICE nodes. Q is fluctuated once per sweep, at t = 1.
    steps = function(problem) {
      lapply(problem$tau:1, function(t) {
        list(id = t, t = t, fluctuate_q = (t == 1L))
      })
    },

    psi_from_state = function(problem, state, detach = TRUE) {
      psi <- state$mu$nodes[[1]]
      if(detach) {
        psi <- psi$detach()$clone()
        psi$requires_grad_(TRUE)
      }
      psi
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      k <- problem$K
      p <- problem$p
      HA_t <- problem$aux$HA_node[, , step$t]
      NablaLdot <- torch::torch_stack(
        batched_NablaLdot(problem$Lm_fn, psi, state$beta,
                          problem$design_matrix, p, k),
        dim = 1
      )
      d_all <- torch::torch_matmul(Minv, NablaLdot$permute(c(2, 1, 3))) # (n, p, k)
      HA_perm <- HA_t$t()$reshape(c(n, k, 1)) # (n, k, 1)

      list(arms = d_all$transpose(2, 3) * HA_perm)
    },

    fluctuation_objective = function(epsilon, problem, step, state, clever, K_Q) {
      mu_node <- state$mu$nodes[[step$t]]
      mu_next <- state$mu$nodes[[step$t + 1]]
      target <- 0
      for(j in seq_len(problem$K)) {
        pred_j <- if(tmle_linear) {
          mu_node[, j] + clever$arms[, j, ]$matmul(epsilon)
        }
        else {
          mu_node[, j]$logit() + clever$arms[, j, ]$matmul(epsilon)
        }
        target <- target + tmle_loss(pred_j, mu_next[, j])
      }
      Qf <- Q_fluctuation(epsilon, K_Q, state$Q)
      target - state$Q$log()$sum()
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q) {
      t <- step$t
      update <- state$mu$nodes[[t]]$detach()$clone()
      for(j in seq_len(problem$K)) {
        update[, j] <- if(tmle_linear) {
          state$mu$nodes[[t]][, j] + clever$arms[, j, ]$matmul(epsilon)
        }
        else {
          torch::torch_sigmoid(
            state$mu$nodes[[t]][, j]$logit() + clever$arms[, j, ]$matmul(epsilon)
          )
        }
      }
      state$mu$nodes[[t]] <- update$detach()
      if(isTRUE(step$fluctuate_Q)) {
        state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
      }
      state$initial <- FALSE
      state
    },

    # Telescoping ICE residual:
    #  Delta[i, j] = sum_t W[j, i, t] (mu_{t+1} - mu_t)[j, i] / g_{0:t}[j, i]
    delta = function(problem, state) {
      n <- problem$n
      k <- problem$K
      tau <- problem$tau
      nodes <- lapply(state$mu$nodes, as.matrix)
      resid <- array(0, dim = c(k, n, tau))
      for(t in seq_len(tau)) resid[, , t] <- t(nodes[[t + 1L]] - nodes[[t]])
      resid <- resid / problem$aux$pi_cumprod * problem$aux$W
      stopifnot(all(is.finite(resid)))
      as_float_tensor(t(apply(resid, c(1, 2), sum))) # (n, k)
    },

    extra_result = function(problem, state) {
      list(tau = problem$tau, regimes = problem$aux$regimes)
    }
  )
}
