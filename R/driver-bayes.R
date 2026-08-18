#' Initial proposal covariance for the generalized posterior
#'
#' The Bayesian submodel is set up so the score at eps = 0 equals
#' the efficient influence function (condition B1), so that
#' dbeta/deps = P_0\[lambda*(P_0)\] = Sigma_D at convergence.
#' Combined with the BvM, which states Cov(beta | O) ~ Sigma_D / n, this gives
#'
#'   Cov(eps | O) ~ Sigma_D^{-1} / n = solve(crossprod(eif)),
#'
#' because \hat{Sigma_D} = crossprod(eif) / n. Optimal scaling
#' for random-walk Metropolis suggests the proposal covariance to be
#' (ell^2 / p) * Cov(target), with ell around 2.38.
#'
#' @param eif n x p matrix, the TMLE EIF
#' @return p x p proposal covariance matrix
#' @noRd
bayes_proposal_scale <- function(eif, ell = 2.38) {
  eif <- as.matrix(eif)
  p <- ncol(eif)

  A <- crossprod(eif)
  S <- tryCatch(solve(A), error = function(e) NULL)

  if(is.null(S) || !all(is.finite(S)) || min(diag(S)) <= 0) {
    warning("Could not derive proposal scale from the influence function.",
            "Falling back to a diagonal proposal; consider supplying `scale` explicitly.", call. = FALSE)

    d <- diag(A)
    d[d <= 0] <- NA_real_
    S <- diag(1 / d, nrow = p)
    if(anyNA(S)) return(diag(1e-3, nrow = p))
  }

  # Symmetrise
  S <- (S + t(S)) / 2
  (ell^2 / p) * S
}

#' Scales clever covariates for Bayesian log likelihood
#' @noRd
scale_bayes_clever <- function(problem, spec, state, clever) {
  clever_scale <- spec$bayes_clever_scale(problem, state)
  if(is.null(clever_scale)) return(clever)
  n <- problem$n
  K <- problem$K
  for(name in names(clever_scale)) {
    cc <- clever[[name]]
    sc <- if(cc$dim() == 2L) {
      clever_scale[[name]]$reshape(c(n, 1))
    }
    else {
      clever_scale[[name]]$reshape(c(n, K, 1))
    }
    clever[[name]] <- cc * sc
  }
  clever
}

#' Jacobian of eps -> beta under the Bayes fluctuation model
#'
#' J(eps) = sum_k dB_dpsi\[k\]' %*% (dpsi_k/deps) + dB_dQ' %*% (dQ/deps)
#'
#' @noRd
bayes_jacobian <- function(problem, psi, Q_base, Q_eps, beta, dpsi_list, epsilon, K_Q, include_Q = TRUE) {
  p <- problem$p
  J_psi <- dB_dpsi(problem$Lm_fn, psi, beta, problem$design_matrix, Q_eps, p) # (K, n, p)
  jac <- torch::torch_zeros(c(p, p))
  for(k in seq_len(problem$K)) {
    jac <- jac + J_psi[k, , ]$transpose(1, 2)$matmul(dpsi_list[[k]])
  }
  if(!include_Q) return(jac)

  Qc <- Q_eps$detach()$clone()$requires_grad_(TRUE)
  J_Q <- dB_dQ(problem$Lm_fn, psi, beta, problem$design_matrix, Qc, p)
  jac <- jac + J_Q$transpose(1, 2)$matmul(dQ_fluctuation_depsilon(epsilon, K_Q, Q_base))
}

#' Generalized Bayesian log-density of the fluctuation parameter
#'
#' The generalized posterior is
#'   log pi(eps | O) = log f(O | eps) + log pi_beta(B(eps)) + log|det J(eps)|
#'
#' Returns the list expected by AdaptMCMC: the first element is the
#' log-density, remaining elements are recorded in `extra.values`
#' @noRd
bayes_log_density <- function(epsilon, problem, spec, state, clever, K_Q, condvar, prior) {
  step <- spec$steps(problem)[[1]]

  clever_b <- scale_bayes_clever(problem, spec, state, clever)

  st <- spec$apply_update(problem, step, state, epsilon, clever_b, K_Q)
  psi <- spec$psi_from_state(problem, st)
  beta <- B(problem$Lm_fn, psi, problem$design_matrix, st$Q, problem$p)

  loglik <- as.numeric(spec$bayes_loglik(epsilon, problem, state, clever_b, st$Q, condvar))
  dpsi <- spec$dpsi_depsilon(problem, st, clever, epsilon)
  jac <- bayes_jacobian(problem, psi, state$Q, st$Q, beta, dpsi, epsilon, K_Q)

  target <- loglik + prior(as.numeric(beta)) + as.numeric(jac$det()$abs()$log())

  list(log.density = as.numeric(target), beta = as.numeric(beta))
}

#' @noRd
run_bayes_tmle <- function(problem, spec, fit, control, eif = NULL) {
  if(!spec$supports_bayes) {
    stop("Generalized Bayesian TMLE is not supported for estimand '", spec$estimand, "'.", call. = FALSE)
  }
  if(!fit$converged) {
    warning("TMLE did not converge; skipping the generalized Bayesian step.")
    return(list(samples = NULL, acc_rate = NA_real_))
  }

  p <- problem$p
  final <- fit$final
  condvar <- if(is.null(problem$nuisance_estimates$condvar)) NULL else as_float_tensor(problem$nuisance_estimates$condvar)

  log_dens <- function(epsilon) {
    bayes_log_density(torch::torch_tensor(epsilon), problem, spec, fit$state, final$clever, final$K_Q, condvar, control$prior)
  }

  total <- control$warmup + control$draws
  keep <- (control$warmup + 1L):total # Discard warmup draws

  beta_samples <- array(NA_real_, dim = c(control$chains, control$draws, p))
  epsilon_samples <- array(NA_real_, dim = c(control$chains, control$draws, p + 1L))
  acc <- numeric(control$chains)

  scale <- control$scale
  if(is.null(scale)) {
    if(is.null(eif)) {
      stop("Automatic proposal scaling requires the TMLE influence function.", call. = FALSE)
    }
    scale <- bayes_proposal_scale(eif)
  }
  else if(!is.matrix(scale)) {
    scale <- diag(rep_len(scale, p), nrow = p)
  }

  run_chain <- function() {
    adaptMCMC::MCMC(
      log_dens, n = total, init = as.numeric(final$epsilon),
      adapt = control$warmup, acc.rate = control$acc_rate, scale = scale
    )
  }

  for(chain in seq_len(control$chains)) {
    mcmc <- if(is.null(control$seed)) {
      run_chain()
    }
    else {
      withr::with_seed(control$seed + chain, run_chain())
    }

    beta_all <- matrix(unlist(mcmc$extra.values), ncol = p, nrow = total, byrow = TRUE)
    beta_samples[chain, , ] <- beta_all[keep, , drop = FALSE]
    epsilon_samples[chain, , ] <- cbind(mcmc$samples, mcmc$log.p)[keep, , drop = FALSE]
    acc[chain] <- mcmc$acceptance.rate
  }

  vars <- beta_variable_names(problem, control$labels)

  list(
    draws = new_draws_array(beta_samples, vars, terms = problem$terms),
    epsilon_draws = new_draws_array(
      epsilon_samples, c(sprintf("epsilon[%d]", seq_len(p)), "lp__")
    ),
    diagnostics = bayes_diagnostics(new_draws_array(beta_samples, vars)),
    acc_rate = mean(acc),
    acc_rate_by_chain = acc,
    scale = scale,
    scale_auto = is.null(control$scale),
    warmup = control$warmup
  )
}
