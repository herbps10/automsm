
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

#' Resolve the compact support of the epsilon submodel
#'
#' @noRd
resolve_bayes_support <- function(problem, control, K_Q, epsilon_init) {
  p <- problem$p
  kmax <- max(abs(as.numeric(K_Q)))

  eps_max <- control$eps_max
  if(is.null(eps_max)) {
    eps_max <- if(kmax > 0) 10 / (2 * p * kmax) else Inf
  }

  min_ess <- control$min_ess %||% (2 * p)

  init_max <- max(abs(as.numeric(epsilon_init)))
  if(is.finite(eps_max) && init_max >= eps_max) {
    warning("The converged TMLE epsilon (max|eps| = ", signif(init_max, 3),
            ") lies outside the resolved eps_max = ", signif(eps_max, 3),
            ". Widening the box so the sampler starts in support.", call. = FALSE)
    eps_max <- 2 * init_max
  }
  list(eps_max = eps_max, min_ess = min_ess, K_Q_max = kmax)
}

#' Scales clever covariates for Bayesian log likelihood
#' @noRd
scale_bayes_clever <- function(problem, spec, state, clever, linear) {
  clever_scale <- spec$bayes_clever_scale(problem, state, linear)
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
bayes_jacobian <- function(problem, psi, Q_base, Q_eps, beta, dpsi_list, epsilon, K_Q, include_Q = TRUE, tol = 1e-10) {
  p <- problem$p
  H <- dobjective_dbeta(problem$Lm_fn, psi, beta, problem$design_matrix, Q_eps, p)
  check <- invert_objective_hessian(H, tol = tol)

  # Check for singular Hessian
  if(is.null(check)) return(NULL)

  J_psi <- dB_dpsi(problem$Lm_fn, psi, beta, problem$design_matrix, Q_eps, p, Hinv = check$inv) # (K, n, p)
  jac <- torch::torch_zeros(c(p, p))
  for(k in seq_len(problem$K)) {
    jac <- jac + J_psi[k, , ]$transpose(1, 2)$matmul(dpsi_list[[k]])
  }
  if(!include_Q) return(jac)

  Qc <- Q_eps$detach()$clone()$requires_grad_(TRUE)
  J_Q <- dB_dQ(problem$Lm_fn, psi, beta, problem$design_matrix, Qc, p, Hinv = check$inv)
  jac <- jac + J_Q$transpose(1, 2)$matmul(dQ_fluctuation_depsilon(epsilon, K_Q, Q_base))
}

#' @noRd
bayes_jacobian_wls <- function(pre, sol, psi, Q, dpsi_list, dQ_deps = NULL) {
  n <- pre$n
  K <- pre$K
  p <- pre$p
  Ai <- chol2inv(sol$chol)
  Qv <- as.numeric(Q)
  stopifnot(length(Qv) == n, length(dpsi_list) == K)

  w <- rep(Qv, times = K) * pre$hrep
  r <- as.numeric(psi) - as.vector(pre$X %*% sol$beta)

  D <- do.call(rbind, lapply(dpsi_list, function(d) {
    if(inherits(d, "torch_tensor")) as.matrix(d$detach()) else as.matrix(d)
  }))
  stopifnot(identical(dim(D), c(n * K, p)))
  M <- crossprod(pre$X * w, D)

  if(!is.null(dQ_deps)) {
    Xw <- pre$X * (pre$hrep * r)
    G <- matrix(0, n, p)
    for(k in seq_len(K)) G <- G + Xw[(k - 1L) * n + seq_len(n), , drop = FALSE]
    M <- M + crossprod(G, dQ_deps)
  }

  Ai %*% M
}

#' Jacobian of the Q-fluctuation in base R
#' @noRd
dQ_deps_wls <- function(K_mat, Q_eps) {
  Qt <- as.numeric(Q_eps)
  n <- length(Qt)
  cs <- as.vector(crossprod(K_mat, Qt))
  Qt * (K_mat - rep(cs, each = n))
}

#' @noRd
log_abs_det <- function(J) {
  if(inherits(J, "torch_tensor")) return(as.numeric(J$det()$abs()$log()))
  d <- determinant(J, logarithm = TRUE)
  if(!is.finite(d$modulus)) return(NA_real_)
  as.numeric(d$modulus)
}

#' Generalized Bayesian log-density of the fluctuation parameter
#'
#' The generalized posterior is
#'   log pi(eps | O) = log f(O | eps) + log pi_beta(B(eps)) + log|det J(eps)|
#'
#' Returns the list expected by AdaptMCMC: the first element is the
#' log-density, remaining elements are recorded in `extra.values`
#' @noRd
bayes_log_density <- function(epsilon, problem, spec, state, clever, K_Q, condvar, prior, linear, beta_init = NULL, control = NULL, tally = NULL) {
  p <- problem$p

  reject <- function(reason) {
    if(!is.null(tally)) tally$n[[reason]] <- (tally$n[[reason]] %||% 0L) + 1L
    list(log.density = -Inf, beta = rep(NA_real_, p))
  }

  eps <- as.numeric(epsilon)
  if(!all(is.finite(eps))) return(reject("nonfinite_eps"))
  if(max(abs(eps)) > (control$eps_max %||% Inf)) return(reject("eps_max"))

  step <- spec$steps(problem)[[1]]

  clever_b <- scale_bayes_clever(problem, spec, state, clever, linear)

  st <- spec$apply_update(problem, step, state, epsilon, clever_b, K_Q, linear)

  # Effective sample size of Q_eps: n at uniform, 1 at a point mass
  ess <- 1 / sum(as.numeric(st$Q)^2)
  if(ess < (control$min_ess %||% 0)) return(reject("low_ess"))

  psi <- spec$psi_from_state(problem, st)

  fast <- problem$B_wls
  fast_ok <- !is.null(fast) && !identical(getOption("automsm.B_wls"), FALSE)
  psi_mat <- as.matrix(psi$detach())
  if(fast_ok) {
    fast_sol <- B_wls(fast, as.matrix(psi), as.numeric(st$Q))
    if(is.null(fast_sol)) return(reject("singular_hessian"))
    beta_num <- fast_sol$beta
  }
  else {
    beta <- B(problem$Lm_fn, psi, problem$design_matrix, st$Q, problem$p, init = beta_init)
    beta_num <- as.numeric(beta)
  }

  if(!all(is.finite(beta_num))) return(reject("nonfinite_beta"))

  loglik <- as.numeric(spec$bayes_loglik(epsilon, problem, state, clever_b, st$Q, condvar, linear = linear))
  dpsi <- spec$dpsi_depsilon(problem, st, clever_b, epsilon, linear = linear)

  if(fast_ok) {
    dQ <- dQ_deps_wls(fast$K_Q, st$Q)
    jac <- bayes_jacobian_wls(fast, fast_sol, as.matrix(psi), st$Q, dpsi, dQ)
  }
  else {
    jac <- bayes_jacobian(problem, psi, state$Q, st$Q, beta, dpsi, epsilon, K_Q)
  }

  if(is.null(jac)) return(reject("singular_hessian"))

  ld <- log_abs_det(jac)
  if(!is.finite(ld)) return(reject("singular_jacobian"))

  target <- loglik + prior(beta_num) + ld
  if(!is.finite(target)) return(reject("nonfinite_target"))

  list(log.density = as.numeric(target), beta = beta_num)
}

#' @noRd
run_bayes_tmle <- function(problem, spec, fit, control, eif = NULL) {
  if(!spec$supports_bayes) {
    stop("Generalized Bayesian TMLE is not supported for estimand '", spec$estimand, "'.", call. = FALSE)
  }
  if(!fit$converged) {
    warning("TMLE did not converge; skipping the generalized Bayesian step.")
    return(NULL)
  }

  p <- problem$p
  final <- fit$final
  condvar <- if(is.null(problem$nuisance_estimates$condvar)) NULL else as_float_tensor(problem$nuisance_estimates$condvar)

  checkmate::assert_number(control$prior(rep(0, p)), finite = TRUE, .var.name = "bayes_control(prior)")

  support <- resolve_bayes_support(problem, control, final$K_Q, final$epsilon)
  control <- utils::modifyList(control, support[c("eps_max", "min_ess")])

  beta_init <- as.numeric(fit$state$beta)


  if(!is.null(problem$B_wls)) problem$B_wls$K_Q <- as.matrix(fit$final$K_Q)

  tally <- new_bayes_tally()
  eval_density <- function(epsilon, tally_env) {
    bayes_log_density(
      torch::torch_tensor(epsilon), problem, spec, fit$state,
      final$clever, final$K_Q, condvar,
      prior = control$prior, linear = control$linear, control = control, tally = tally_env, beta_init = beta_init
    )
  }
  log_dens <- function(epsilon) eval_density(epsilon, tally)


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

  probe_tally <- new_bayes_tally()
  ld0 <- eval_density(final$epsilon, probe_tally)$log.density
  if(!is.finite(ld0)) {
    stop("Generalizd posterior log-density is not finite at the converged TMLE epsilon.", call. = FALSE)
  }

  total <- control$warmup + control$draws
  keep <- (control$warmup + 1L):total # Discard warmup draws

  beta_samples <- array(NA_real_, dim = c(control$chains, control$draws, p))
  epsilon_samples <- array(NA_real_, dim = c(control$chains, control$draws, p + 1L))
  acc <- numeric(control$chains)
  rej_by_chain <- vector("list", control$chains)


  run_chain <- function() {
    adaptMCMC::MCMC(
      log_dens, n = total, init = as.numeric(final$epsilon),
      adapt = control$warmup, acc.rate = control$acc_rate, scale = scale
    )
  }

  for(chain in seq_len(control$chains)) {
    before <- tally_counts(tally)

    mcmc <- if(is.null(control$seed)) {
      run_chain()
    }
    else {
      withr::with_seed(control$seed + chain, run_chain())
    }

    rej_by_chain[[chain]] <- tally_diff(tally_counts(tally), before)

    beta_all <- matrix(unlist(mcmc$extra.values), ncol = p, nrow = total, byrow = TRUE)
    beta_samples[chain, , ] <- beta_all[keep, , drop = FALSE]
    epsilon_samples[chain, , ] <- cbind(mcmc$samples, mcmc$log.p)[keep, , drop = FALSE]
    acc[chain] <- mcmc$acceptance.rate
  }

  rejected <- tally_counts(tally)
  if(sum(rejected) > 0.5 * control$chains * total) {
    warning("More than half of all proposals fell outside the epsilon support ",
            paste(names(rejected), rejected, sep = "=", collapse = ", "),
            "). Tune `scale`/`warmup.", call. = FALSE)
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
    rejected = rejected,
    rejected_by_chain = rej_by_chain,
    scale = scale,
    scale_auto = is.null(control$scale),
    eps_max = control$eps_max, min_ess = control$min_ess, K_Q_max = support$K_Q_max,
    warmup = control$warmup,
    linear = control$linear,
    fluctuation = control$fluctuation
  )
}

#' Mutable counter for rejected MCMC proposals
#' @noRd

new_bayes_tally <- function() {
  t <- new.env(parent = emptyenv())
  t$n <- list()
  t
}

tally_counts <- function(tally) {
  if(is.null(tally) || length(tally$n) == 0L) return(integer(0))
  unlist(tally$n)
}

tally_get <- function(counts, name) {
  out <- integer(length(name))
  names(out) <- name
  common <- intersect(name, names(counts))
  out[common] <- as.integer(counts[common])
  out
}

tally_diff <- function(after, before) {
  name <- union(names(after), names(before))
  if(length(name) == 0L) return(integer(0))
  tally_get(after, name) - tally_get(before, name)
}
