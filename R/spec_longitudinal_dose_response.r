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

    psi_from_state = function(problem, state) {
      psi <- state$mu$nodes[[1]]
      psi <- psi$detach()$clone()
      psi$requires_grad_(TRUE)
      psi
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      K <- problem$K
      d <- clever_directions(problem, psi, state$beta, Minv)
      list(arms = d * problem$aux$HA_node[, , step$t]$t()$reshape(c(n, K, 1)))
    },

    mu_loss = function(epsilon, problem, step, state, clever) {
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
      target
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q) {
      t <- step$t
      lin <- clever$arms$matmul(epsilon)
      node <- state$mu$nodes[[t]]
      state$mu$nodes[[t]] <- if(tmle_linear) {
        (node + lin)$detach()
      }
      else {
        torch::torch_sigmoid(node$logit() + lin)$detach()
      }
      if(isTRUE(step$fluctuate_Q)) state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
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

    nuisance_contract = function(problem, bayes_enabled) {
      n <- problem$n
      K <- problem$K
      tau <- problem$tau
      b <- outcome_bounds(problem)
      list(
        fields = list(
          nuisance_field("pi", c(n, tau), lower = 0, upper = 1),
          nuisance_field("mu", c(K, n, tau + 1), lower = b$lo, upper = b$hi, severity = "warning")
        ),
        checks = list(
          function(nu, problem) {
            Y <- as.numeric(problem$Yt)
            ok <- vapply(seq_len(K), function(j) {
              max(abs(nu$mu[j, , problem$tau + 1L] - Y)) < 1e-6
            }, logical(1))
            if(all(ok)) TRUE else "nuisance$mu[, , tau + 1] must equal Y for every regime."
          }
        )
      )
    }
  )
}
