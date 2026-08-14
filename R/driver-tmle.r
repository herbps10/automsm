run_tmle <- function(problem, spec, control, state) {
  steps <- spec$steps(problem)
  tol <- control$tol %||% spec$tol
  converged <- TRUE
  iter <- 0L
  cache <- NULL

  for(it in seq_len(control$maxiter)) {
    iter <- it
    eps_list <- vector("list", length(steps))

    for(s in seq_along(steps)) {
      step <- steps[[s]]

      psi <- spec$psi_from_state(problem, state, detach = TRUE)
      Minv <- normalizing_matrix(problem$Lm_fn, psi, state$beta,
                                 problem$design_matrix, state$Q, problem$p)
      K_Q <- calculate_K(problem$Lm_fn, psi, state$beta,
                         problem$design_matrix, Minv)
      clever <- spec$make_clever(problem, step, state, psi, Minv)
      eps <- tmle_mle(problem$p, function(epsilon) {
        spec$fluctuation_objective(epsilon, problem, step, state, clever, K_Q)
      })

      if(spec$nan_guard && any(is.nan(as.numeric(eps)))) {
        warning("TMLE failed to converge.")
        converged <- FALSE
        break
      }

      eps_list[[s]] <- eps

      cache <- list(step = step, psi = psi, Minv = Minv, K_Q = K_Q,
                    clever = clever, epsilon = eps)

      state <- spec$apply_update(problem, step, state, eps, clever, K_Q)

    }
    if(!converged) break

    state <- update_beta(problem, spec, state)

    if(abs(spec$sweep_criterion(eps_list)) < tol) break
  }

  list(state = state, converged = converged, iter = iter,
       cache = cache, steps = steps)
}

update_beta <- function(problem, spec, state) {
  psi <- spec$psi_from_state(problem, state, detach = TRUE)
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

  psi <- spec$psi_from_state(problem, fit$state, detach = TRUE)

  est <- B(problem$Lm_fn, psi, problem$design_matrix, fit$state$Q, problem$p)
  Delta <- spec$delta(problem, fit$state)
  eifm <- eif(problem$Lm_fn, psi, est, problem$design_matrix, fit$state$Q, Delta, problem$p)
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

run_bayes_tmle_shim <- function(problem, spec, fit, tmle_linear, draws, chains, prior) {
  p <- problem$p
  state <- fit$state
  cache <- fit$cache
  samples <- array(dim = (c(chains, draws, p)))
  acc <- 0

  log_dens <- function(epsilon) {
    old_cate_fluctuation_model(
      torch::torch_tensor(epsilon),
      problem$Lm_fn,
      state$mu$obs, state$mu$a0, state$mu$a1,
      cache$clever$obs, cache$clever$a0, cache$clever$a1,
      cache$K_Q, state$Q, problem$Yt, problem$design_matrix,
      condvar = torch::torch_tensor(problem$nuisance$condvar),
      bayes = TRUE,
      bayes_prior = prior,
      tmle_loss = spec$tmle_loss
    )
  }

  for(chain in seq_len(chains)) {
    mcmc <- adaptMCMC::MCMC(log_dens, n = draws, init = as.numeric(cache$epsilon),
                            adapt = TRUE, acc.rate = 0.3, scale = rep(1e-3, p))

    samples[chain, , ] <- matrix(unlist(mcmc$extra.values), ncol = p, nrow = draws, byrow = TRUE)

    acc <- acc + mcmc$acceptance.rate / chains
  }
  list(samples = samples, acc_rate = acc)
}
