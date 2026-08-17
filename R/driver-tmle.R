run_tmle <- function(problem, spec, control, state) {
  steps <- spec$steps(problem)
  tol <- control$tol %||% spec$tol
  converged <- TRUE
  iter <- 0L

  for(it in seq_len(control$maxiter)) {
    iter <- it
    eps_list <- vector("list", length(steps))

    for(s in seq_along(steps)) {
      step <- steps[[s]]

      psi <- spec$psi_from_state(problem, state)
      Minv <- normalizing_matrix(problem$Lm_fn, psi, state$beta,
                                 problem$design_matrix, state$Q, problem$p)
      K_Q <- calculate_K(problem$Lm_fn, psi, state$beta,
                         problem$design_matrix, Minv, problem$batched_beta)
      clever <- spec$make_clever(problem, step, state, psi, Minv)

      obj <- function(epsilon) {
        target <- spec$mu_loss(epsilon, problem, step, state, clever)
        if(isTRUE(steps$fluctuate_Q)) {
          target <- target - Q_fluctuation(epsilon, K_Q, state$Q)$log()$sum()
        }
        target
      }

      eps <- tmle_mle(problem$p, obj)

      if(spec$nan_guard && any(is.nan(as.numeric(eps)))) {
        warning("TMLE failed to converge.")
        converged <- FALSE
        break
      }

      eps_list[[s]] <- eps

      state <- spec$apply_update(problem, step, state, eps, clever, K_Q)

    }
    if(!converged) break

    state <- update_beta(problem, spec, state)

    if(abs(spec$sweep_criterion(eps_list)) < tol) break
  }

  if(converged) {
    psi_f <- spec$psi_from_state(problem, state)
    Minv_f <- normalizing_matrix(problem$Lm_fn, psi_f, state$beta,
                                 problem$design_matrix, state$Q, problem$p)

    final <- list(
      psi = psi_f, Minv = Minv_f,
      K_Q = calculate_K(problem$Lm_fn, psi_f, state$beta,
                        problem$design_matrix, Minv_f, problem$batched_beta),
      clever = spec$make_clever(problem, steps[[length(steps)]], state, psi_f, Minv_f),
      epsilon = eps_list[[length(steps)]]
    )
  } else {
    final <- NULL
  }

  list(state = state, converged = converged, iter = iter, steps = steps, final = final)
}

update_beta <- function(problem, spec, state) {
  psi <- spec$psi_from_state(problem, state)
  beta <- B(problem$Lm_fn, psi, problem$design_matrix, state$Q, problem$p)$detach()$clone()
  beta$requires_grad_(TRUE)
  state$beta <- beta
  state
}

finalize_tmle <- function(problem, spec, fit) {
  n <- problem$n
  p <- problem$p

  if(!fit$converged) {
    return(list(
      est = rep(NA_real_, p),
      se = rep(NA_real_, p),
      lower = rep(NA_real_, p),
      upper = rep(NA_real_, p),
      eif = matrix(NA_real_, nrow = n, ncol = p),
      converged = FALSE,
      iter = fit$iter
    ))
  }

  psi <- spec$psi_from_state(problem, fit$state)

  est <- B(problem$Lm_fn, psi, problem$design_matrix, fit$state$Q, problem$p)
  Delta <- spec$delta(problem, fit$state)
  eifm <- eif(problem$Lm_fn, psi, est, problem$design_matrix, fit$state$Q, Delta, problem$p, problem$batched_beta)
  se <- apply(eifm, 2, stats::sd) / sqrt(n)

  list(
    est = as.numeric(est),
    se = as.numeric(se),
    lower = as.numeric(est + stats::qnorm(0.025) * se),
    upper = as.numeric(est + stats::qnorm(0.975) * se),
    eif = eifm,
    converged = TRUE,
    iter = fit$iter
  )
}
