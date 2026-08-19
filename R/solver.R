# Fluctuation solvers
#
# For each estimand and fluctuation type (linear or logistic), the objective is a GLM in epsilon:
#
#   eta = offset + X eps        X = clever covariates
#   f(eps) = sum_i dev(eta_i, y_i)     [+ q_penalty((eps) if step$fluctuate_Q]
#
#   linear = TRUE: dev = (eta - y)^2               (Gaussian, identity link)
#   linear = FALSE: dev = log(1 + e^eta) - y eta   (Bernoulli, logit link)
#
# This is a classical TMLE update. Gradient and Hessian are available in closed form,
# and no autodiff is needed.

#' @noRd
fluctuation_family <- function(linear) {
  if(linear) {
    list(
      name = "gaussian",
      dev = function(eta, y) 0.5 * (eta - y)^2,
      grad = function(eta, y) eta - y,
      hess_w = function(eta, y) rep_len(1, length(eta))
    )
  }
  else {
    list(
      name = "binomial",
      dev = function(eta, y) ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))) - y * eta,
      grad = function(eta, y) stats::plogis(eta) - y,
      hess_w = function(eta, y) {
        m <- stats::plogis(eta)
        m * (1 - m)
      }
    )
  }
}

#' Analytic gradient and Hessian of the -log Q_eps penalty
#' @noRd
q_penalty <- function(K_Q, Q, n, check_clamp = TRUE) {
  shift <- function(eps) {
    u <- as.vector(K_Q %*% eps)
    u - max(u)
  }

  list(
    value = function(eps) {
      u <- shift(eps)
      if(check_clamp && min(u) < -60) {
        warning("Q-fluctuation underflow clamp is active; the analytic penalty no longer matches Q_fluctuation() exactly.", call. = FALSE)
      }
      -sum(u) - sum(log(Q)) + n * log(sum(exp(u) * Q))
    },
    grad_hess = function(eps) {
      u <- shift(eps)
      w <- exp(u) * Q
      Qt <- w / sum(w)
      KQt <- as.vector(crossprod(K_Q, Qt))
      list(g = -colSums(K_Q) + n * KQt,
           H = n * (crossprod(K_Q * Qt, K_Q) - outer(KQt, KQt)))
    }
  )
}

#' Objective value of a weighted GLM fluctuation
#' @noRd
glm_objective <- function(glm, family, eps, qpen = NULL) {
  w <- glm$weights %||% 1
  eta <- as.vector(glm$offset + glm$X %*% eps)
  v <- sum(w * family$dev(eta, glm$target))
  if(!is.null(qpen)) v <- v + qpen$value(eps)
  if(!is.finite(v)) return(sqrt(.Machine$double.xmax))
  v
}

#' Analytic gradient of glm_objective()
#' @noRd
glm_gradient <- function(glm, family, eps, qpen = NULL) {
  w <- glm$weights %||% 1
  eta <- as.vector(glm$offset + glm$X %*% eps)
  g <- as.vector(crossprod(glm$X, w * family$grad(eta, glm$target)))
  if(!is.null(qpen)) g <- g + qpen$grad_hess(eps)$g
  g
}

#' Gradient-only L-BFGS-B fallback
#'
#' Only uses the graident. tmle_newton() fails when when Hessian collapses,
#' so a fallback that also uses the inverted Hessian would also fail.
#'
#' @noRd
tmle_lbfgsb <- function(p, glm, family, qpen = NULL, eps_max = Inf,
                        maxit = 200L, factr = 1e2) {
  fn <- function(e) glm_objective(glm, family, e, qpen)
  gr <- function(e) glm_gradient(glm, family, e, qpen)
  bound <- if(is.finite(eps_max)) eps_max else Inf

  out <- stats::optim(rep(0, p), fn = fn, gr = gr, method = "L-BFGS-B",
                      lower = rep(-bound, p), upper = rep(bound, p),
                      control = list(maxit = maxit, factr = factr))

  list(epsilon = out$par,
       iter = as.integer(out$counts[["function"]]),
       grad_inf = max(abs(gr(out$par))),
       value = out$value,
       truncated = is.finite(eps_max) && max(abs(out$par)) >= eps_max * (1 - 1e-8),
       singular = FALSE,
       convergence = out$convergence,
       message = out$message %||% "")
}

#' Solve a weighted GLM fluctuation step
#'
#' Newton/Iteratively Reweighted Least Squares with L-BFGS fallback
#' @noRd
solve_fluctuation_glm <- function(p, glm, family, control, label = "") {
  f0 <- sum((glm$weights %||% 1) * family$dev(as.vector(glm$offset), glm$target))

  acceptable <- function(f) {
    is.finite(f) && f <= f0 + control$obj_tol * (abs(f0) + 1L)
  }

  try_solver <- function(f, name) {
    out <- tryCatch(
      f(),
      error = function(e) {
        warning("Solver '", name, "' failed at ", label, " (", conditionMessage(e), "); falling back to L-BFGS.", call. = FALSE)
        NULL
      }
    )
    if(is.null(out)) return(NULL)
    out$solver <- name
    if(!acceptable(out$value)) {
      warning("Solver '", name, "' did not decrease the objective at ", label,
              " (f(0) = ", signif(f0, 6), ", f(eps) = ", signif(out$value, 6), ").", call. = FALSE)
      return(NULL)
    }
    out
  }

  out <- try_solver(function() tmle_newton(p, glm, family, qpen = NULL, eps_max = control$eps_max, ridge = control$ridge), "newton")

  if(!is.null(out) && isTRUE(out$singular)) {
    warning("Fluctuation Hessian collapsed at ", label,
            "; the submodel may be separated or the design degenerate. ",
            "Retrying with the gradient-only solver.", call. = FALSE)
    out <- NULL
  }

  if(is.null(out)) {
    out <- try_solver(function() tmle_lbfgsb(p, glm, family, qpen = NULL, eps_max = control$eps_max), "lbfgsb")
  }

  # Objective is convex, so eps = 0 must not beat the returned iterate
  if(is.null(out)) {
    warning("No solver decreased the objective at ", label, "; taking eps = 0 for this node.", call. = FALSE)
    out <- list(
      epsilon = rep(0, p), solver = "none", iter = NA_integer_,
      grad_inf = NA_real_, value = f0, truncated = FALSE,
      singular = FALSE, convergence = NA_integer_, message = ""
    )
  }

  out$improved <- acceptable(out$value)
  out$log <- out[c("solver", "iter", "grad_inf", "value", "improved", "truncated", "singular", "convergence")]
  out$log$max_abs_eps <- max(abs(out$epsilon))
  out
}

#' Newton/Iteratively Reweighted Least Squares solver for the fluctuation parameter
#' @noRd
tmle_newton <- function(p, glm, family, qpen = NULL, eps_max = Inf, maxit = 50L, tol = 1e-10, ridge = 1e-10) {
  X <- glm$X
  offset <- glm$offset
  y <- glm$target
  stopifnot(nrow(X) == length(offset), length(offset) == length(y), ncol(X) == p)

  obj <- function(e) {
    v <- sum((glm$weights %||% 1) * family$dev(as.vector(offset + X %*% e), y))
    if(!is.null(qpen)) v <- v + qpen$value(e)
    v
  }

  eps <- rep(0, p)
  f0 <- obj(eps)
  g <- rep(NA_real_, p)
  it <- 0L
  truncated <- FALSE
  singular <- FALSE
  convergence <- FALSE

  for(it in seq_len(maxit)) {
    eta <- as.vector(offset + X %*% eps)
    w <- glm$weights %||% 1
    W <- family$hess_w(eta, y)
    g <- as.vector(crossprod(X, w * family$grad(eta, y)))
    H <- crossprod(X * w * W, X)
    if(!is.null(qpen)) {
      q <- qpen$grad_hess(eps)
      g <- g + q$g
      H <- H + q$H
    }
    H <- H + diag(ridge * max(1, max(abs(H))), p)

    d <- tryCatch(as.vector(solve(H, g)), error = function(e) NULL)
    if(is.null(d) || !all(is.finite(d))) {
      # Collapsed Hessian
      # the submodel is separated or the design is degenerate
      singular <- TRUE
      break
    }

    # Backtracking line search
    tt <- 1
    slope <- sum(g * d)
    cand <- eps
    repeat {
      cand <- eps - tt * d
      f1 <- obj(cand)
      if(is.finite(f1) && f1 <= f0 - 1e-4 * tt * slope) break
      tt <- tt / 2
      if(tt < 1e-12) {
        cand <- eps
        f1 <- f0
        break
      }
    }

    if(max(abs(cand)) > eps_max) {
      cand <- cand * eps_max / max(abs(cand))
      f1 <- obj(cand)
      truncated <- TRUE
    }

    delta <- max(abs(cand - eps))
    eps <- cand
    f0 <- f1
    if(delta < tol || max(abs(g)) < tol) {
      convervence <- TRUE
      break
    }
  }

  list(epsilon = eps, iter = it, grad_inf = max(abs(g)), value = f0,
       truncated = truncated, singular = singular, convergence = convergence )
}

#' L-BFGS fallback
#' @noRd
tmle_mle <- function(p, obj_fn) {
  epsilon <- torch::torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(epsilon, max_iter = 20, line_search_fn = "strong_wolfe")

  for (iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      target <- obj_fn(epsilon)
      target$backward(retain_graph = TRUE)
      target
    })
  }
  optimizer$zero_grad()
  epsilon
}

#' Solve one fluctuation step, with fallback
#' @noRd
solve_fluctuation <- function(problem, spec, control, step, state, clever, K_Q) {
  p <- problem$p

  obj_torch <- function(epsilon) {
    target <- spec$mu_loss(epsilon, problem, step, state, clever, control$linear)

    if(isTRUE(step$fluctuate_Q)) {
      target <- target - Q_fluctuation(epsilon, K_Q, state$Q)$log()$sum()
    }

    target
  }

  value <- function(e) as.numeric(obj_torch(torch::torch_tensor(e)))
  f0 <- value(rep(0, p))
  acceptable <- function(f) {
    is.finite(f) && f <= f0 + control$obj_tol * (abs(f0) + 1L)
  }

  out <- NULL

  if(identical(control$solver, "newton")) {
    g <- spec$fluctuation_glm(problem, step, state, clever)

    gp <- if(isTRUE(step$fluctuate_Q)) {
      q_penalty(as.matrix(K_Q), as.numeric(state$Q), problem$n)
    }
    else {
      NULL
    }

    out <- tryCatch(
      tmle_newton(p, g, fluctuation_family(control$linear), gp, eps_max = control$eps_max, ridge = control$ridge),
      error = function(e) {
        warning("Newton fluctuation solver failed (", conditionMessage(e), "); falling back to L-BFGS.", call. = FALSE)
        NULL
      }
    )

    if(!is.null(out)) {
      out$value_torch <- value(out$epsilon)
      out$solver <- "newton"
      if(!acceptable(out$value_torch)) {
        warning("Newton step did not decrease the torch objective (f(0) = ", signif(f0, 6), ", f(eps) = ", signif(out$value_torch, 6), "); falling back to L-BFGS. ",
                "If this recurs, check spec$fluctuation_glm() against spec$mu_loss().", call. = FALSE)
        out <- NULL
      }

      if(!is.null(out) && isTRUE(out$singular)) {
        warning("Fluctuation Hessian collapsed at step t = ", step$t %||% step$id, "; the submodel may be separated or the clever-covariate design degenerate.", call. = FALSE)
      }
    }
  }

  if(is.null(out)) {
    eps <- as.numeric(tmle_mle(p, obj_torch))
    if(max(abs(eps)) > control$eps_max) {
      eps <- eps * control$eps_max / max(abs(eps))
      out_trunc <- TRUE
    }
    else {
      out_trunc <- FALSE
    }

    out <- list(
      epsilon = eps, solver = "lbfgs", iter = NA_integer_, grad_inf = NA_real_, value_torch = value(eps), truncated = out_trunc, singular = FALSE
    )
  }

  # Monotinicity assertion
  # Objective is convex in epsilon, so eps = 0 must not beat the returned iterate

  out$improved <- acceptable(out$value_torch)
  if(!out$improved) {
    warning("No fluctuation solver decreased the objective at step t = ", step$t %||% step$id,
            " (f(0) = ", signif(f0, 6), ", f(eps) = ", signif(out$value_torch, 6), "); taking eps = 0 for this step.", call. = FALSE)
    out$epsilon <- rep(0, p)
    out$value_torch <- f0
  }

  out$f0 <- f0
  out$epsilon_t <- torch::torch_tensor(out$epsilon)
  out
}
