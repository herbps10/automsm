msm_spec_cate <- function(tmle_linear = TRUE, bayes_linear = TRUE) {
  tmle_loss_linear <- {
    f <- torch::nn_mse_loss(reduction = "sum")
    function(x, y) 0.5 * f(x, y)
  }

  tmle_loss_logistic <- torch::nn_bce_with_logits_loss(reduction = "sum")

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
    #tmle_loss = tmle_loss,

    init_state = function(problem) {
      n <- problem$n
      list(
        mu = list(
          obs = as_float_tensor(problem$nuisance_estimates$mu)$reshape(c(n, 1)),
          a0 = as_float_tensor(problem$nuisance_estimates$mu0)$reshape(c(n, 1)),
          a1 = as_float_tensor(problem$nuisance_estimates$mu1)$reshape(c(n, 1))
        ),
        Q = problem$Q0,
        beta = NULL
      )
    },

    steps = function(problem) list(list(id = 1L, fluctuate_Q = TRUE)),

    psi_from_state = function(problem, state) {
      psi <- state$mu$a1 - state$mu$a0
      psi <- psi$detach()$clone()
      psi$requires_grad_(TRUE)
      psi
    },

    make_clever = function(problem, step, state, psi, Minv) {
      n <- problem$n
      d <- clever_directions(problem, psi, state$beta, Minv)$reshape(c(n, problem$p))
      list(
        obs = d * problem$aux$H$reshape(c(n, 1)),
        a0  = d * problem$aux$H0$reshape(c(n, 1)),
        a1  = d * problem$aux$H1$reshape(c(n, 1))
      )
    },

    fluctuation_offset = function(problem, step, state, linear) {
      m <- state$mu$obs[, 1]
      if(linear) {
        m
      }
      else {
        m$logit()
      }
    },

    fluctuation_glm = function(problem, step, state, clever) {
      list(X = as.matrix(clever$obs),
           offset = as.numeric(clever$offset),
           target = as.numeric(problem$Yt))
    },

    mu_loss = function(epsilon, problem, step, state, clever, linear) {
      p <- clever$offset + clever$obs$matmul(epsilon)
      if(isTRUE(linear)) tmle_loss_linear(p, problem$Yt) else tmle_loss_logistic(p, problem$Yt)
    },

    apply_update = function(problem, step, state, epsilon, clever, K_Q, linear) {
      n <- problem$n

      fluctuate <- function(col, cc) {
        if(isTRUE(linear)) {
          col + cc$matmul(epsilon)
        }
        else {
          clamp_fit(torch::torch_sigmoid(col$logit() + cc$matmul(epsilon)), problem$clamp)
        }
      }

      state$mu <- list(
        obs = fluctuate(state$mu$obs[, 1], clever$obs)$reshape(c(n, 1))$detach(),
        a0 = fluctuate(state$mu$a0[, 1], clever$a0)$reshape(c(n, 1))$detach(),
        a1 = fluctuate(state$mu$a1[, 1], clever$a1)$reshape(c(n, 1))$detach()
      )

      if(isTRUE(step$fluctuate_Q)) {
        state$Q <- Q_fluctuation(epsilon, K_Q, state$Q)
      }

      state
    },

    delta = function(problem, state) {
      (problem$aux$H * (problem$Yt - state$mu$obs[, 1]))$reshape(c(problem$n, 1))
    },

    dpsi_depsilon = function(problem, state, clever, epsilon, linear) {
      n <- problem$n
      if(linear) {
        list(clever$a1 - clever$a0)
      }
      else {
        m1 <- state$mu$a1$reshape(c(n, 1))
        m0 <- state$mu$a0$reshape(c(n, 1))
        list(clever$a1 * m1 * (1 - m1) - clever$a0 * m0 * (1 - m0))
      }
    },

    bayes_clever_scale = function(problem, state, linear) {
      if(linear == FALSE) return(NULL)

      nu <- problem$nuisance_estimates
      list(
        obs = as_float_tensor(nu$condvar),
        a0 = as_float_tensor(nu$condvar0),
        a1 = as_float_tensor(nu$condvar1)
      )
    },

    bayes_loglik = function(epsilon, problem, state, clever, Q_eps, condvar, linear) {
      logQ <- Q_eps$log()$sum()
      if(linear) {
        stopifnot(!is.null(condvar))
        pred <- state$mu$obs[, 1] + clever$obs$matmul(epsilon)
        -((pred - problem$Yt)$pow(2) / (2 * condvar))$sum() - 0.5 * condvar$log()$sum() + logQ
      }
      else {
        lg <- state$mu$obs[, 1]$logit() + clever$obs$matmul(epsilon)
        -tmle_loss_logistic(lg, problem$Yt) + logQ
      }
    },

    nuisance_contract = function(problem, bayes_enabled) {
      n <- problem$n
      b <- outcome_bounds(problem)

      list(
        fields = list(
          nuisance_field("pi", n, lower = 0, upper = 1),
          nuisance_field("mu", n, lower = b$lo, upper = b$hi, severity = "warning"),
          nuisance_field("mu0", n, lower = b$lo, upper = b$hi, severity = "warning"),
          nuisance_field("mu1", n, lower = b$lo, upper = b$hi, severity = "warning"),
          nuisance_field("condvar", n, required = bayes_enabled),
          nuisance_field("condvar0", n, required = bayes_enabled),
          nuisance_field("condvar1", n, required = bayes_enabled)
        ),
        checks = list(
          function(nu, problem) {
            expect <- ifelse(problem$aux$A_obs == 1, nu$mu1, nu$mu0)
            if(max(abs(nu$mu - expect)) < 1e-6) TRUE else "nuisance$mu must equal mu1 where A == 1 and mu0 where A == 0."
          },
          function(nu, problem) {
            expect <- ifelse(problem$aux$A_obs == 1, nu$condvar1, nu$condvar0)
            if(max(abs(nu$condvar - expect)) < 1e-6) TRUE else "nuisance$condvar must equal condvar1 where A == 1 and condvar0 where A == 0."
          }
        )
      )
    }
  )
}
