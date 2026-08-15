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

    psi_from_state = function(problem, state) {
      psi <- state$mu$arms
      psi <- psi$detach()$clone()
      psi$requires_grad_(TRUE)
      psi
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      K <- problem$K
      d <- clever_directions(problem, psi, state$beta, Minv)
      list(
        obs = (d * problem$aux$H$reshape(c(n, K, 1)))$sum(dim = 2),
        arms = d * problem$aux$HA$reshape(c(n, K, 1))
      )
    },

    mu_loss = function(epsilon, problem, step, state, clever) {
      pred <- if(tmle_linear) {
        state$mu$obs + clever$obs$matmul(epsilon)
      }
      else {
        state$mu$obs$logit() + clever$obs$matmul(epsilon)
      }
      tmle_loss(pred, problem$Yt)
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q) {
      obs_lin <- clever$obs$matmul(epsilon)
      arms_lin <- clever$arms$matmul(epsilon)
      if(tmle_linear) {
        obs <- state$mu$obs + obs_lin
        arms <- state$mu$arms + arms_lin
      }
      else {
        obs <- torch::torch_sigmoid(state$mu$obs$logit())
        arms <- torch::torch_sigmoid(state$mu$arms$logit() + arms_lin)
      }
      state$mu <- list(obs = obs$detach(), arms = arms$detach())
      if(isTRUE(step$fluctuate_Q)) state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
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
    },

    nuisance_contract = function(problem, bayes_enabled) {
      n <- problem$n
      K <- problem$K
      b <- outcome_bounds(problem)

      list(
        fields = list(
          nuisance_field("pi", n, lower = 0, upper = 1),
          nuisance_field("mu", n, lower = b$lo, upper = b$hi, severity = "warning"),
          nuisance_field("mu_a", c(n, K), lower = b$lo, upper = b$hi, severity = "warning"),
          nuisance_field("mu_a", c(n, K), lower = 0, 1),
          nuisance_field("condvar", n, required = bayes_enabled && tmle_linear)
        ),
        checks = list(
          function(nu, problem) {
            idx <- cbind(seq_len(problem$n), problem$aux$A_index)
            if(max(abs(nu$pi - as.matrix(nu$pi_a)[idx])) < 1e-6) {
              TRUE
            }
            else {
              paste("nuisance$pi must be the propensity of the observed treatment, i.e. pi[i] == pi_a[i, k] where As[k] == A[i].")
            }
          },
          function(nu, problem) {
            idx <- cbind(seq_len(problem$n), problem$aux$A_index)
            if(max(abs(nu$mu - as.matrix(nu$mu_a)[idx])) < 1e-6) TRUE else "nuisance$mu must equal mu_a[i, k] where As[k] == A[i]."
          }
        )
      )
    }
  )
}
