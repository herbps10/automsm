#' Diagnostics for the EIF estimating equation at current state
#'
#' Returns max_j |Pn\[ D*_j \]| / se_j, which under Theorem 1 is required to be o_P(n^{-1/2}).
#' @noRd
eif_diagnostics <- function(problem, spec, state) {
  psi <- spec$psi_from_state(problem, state)
  est <- state$beta
  Delta <- spec$delta(problem, state)
  eifm <- eif(problem$Lm_fn, psi, est, problem$design_matrix, state$Q, Delta, problem$p, problem$batched_beta)
  se <- apply(eifm, 2, stats::sd) / sqrt(problem$n)
  ratio <- max(abs(colMeans(eifm)) / pmax(se, .Machine$double.eps))
  list(psi = psi, est = est, eif = eifm, se = se, ratio = ratio)
}

#' @noRd
update_beta <- function(problem, spec, state) {
  psi <- spec$psi_from_state(problem, state)
  beta <- B(problem$Lm_fn, psi, problem$design_matrix, state$Q, problem$p)$detach()$clone()
  beta$requires_grad_(TRUE)
  state$beta <- beta
  state
}

#' @noRd
run_tmle <- function(problem, spec, control, state, verbose = FALSE) {
  steps <- spec$steps(problem)
  tol <- control$tol %||% spec$tol
  converged <- TRUE
  iter <- 0L
  status <- "maxiter"
  eps_list <- vector("list", length(steps))
  solver_log <- list()
  diag <- NULL

  for(it in seq_len(control$maxiter)) {
    iter <- it
    sweep_ok <- TRUE

    for(s in seq_along(steps)) {
      step <- steps[[s]]

      psi <- spec$psi_from_state(problem, state)
      Minv <- normalizing_matrix(problem$Lm_fn, psi, state$beta,
                                 problem$design_matrix, state$Q, problem$p)
      K_Q <- calculate_K(problem$Lm_fn, psi, state$beta,
                         problem$design_matrix, Minv, problem$batched_beta)
      clever <- spec$make_clever(problem, step, state, psi, Minv)
      clever$offset <- spec$fluctuation_offset(problem, step, state)

      sol <- solve_fluctuation(problem, spec, control, step, state, clever, K_Q)
      eps <- sol$epsilon_t
      eps_list[[s]] <- eps
      solver_log[[length(solver_log) + 1L]] <- list(
        iter = it, step = s, t = step$t %||% step$id, solver = sol$solver,
        inner = sol$iter, grad_inf = sol$grad_inf, improved = sol$improved,
        truncated = sol$truncated, singular = sol$singular,
        max_abs_eps = max(abs(sol$epsilon))
      )

      if(verbose) {
        cat(sprintf("sweep %2d t=%s solver=%-6s max|eps| = %.3e%s\n",
                    it, step$t %||% step$id, sol$solver, max(abs(sol$epsilon)),
                    if(isTRUE(sol$truncated)) "  [truncated]" else ""))
      }

      if(spec$nan_guard && !all(is.finite(sol$epsilon))) {
        warning("TMLE fluctuation produced non-finite epsilon.", call. = FALSE)
        status <- "nan"
        sweep_ok <- FALSE
        break
      }

      if(isTRUE(sol$singular)) {
        status <- "degenerate"
        sweep_ok <- FALSE
        break
      }

      state <- spec$apply_update(problem, step, state, eps, clever, K_Q)
    }

    if(!sweep_ok) break

    state <- update_beta(problem, spec, state)
    diag <- eif_diagnostics(problem, spec, state)

    done <- if(identical(control$criterion, "eif")) {
      diag$ratio < control$eif_tol
    }
    else {
      abs(spec$sweep_criterion(eps_list)) < tol
    }
    if(verbose) cat(sprintf("  -> solved = %.3e se-units\n", diag$ratio))
    if(done) {
      status <- "converged"
      break
    }
  }

  if(identical(status, "maxiter") && !is.null(diag)) {
    warning("TMLE reached maxiter = ", control$maxiter, " without converging ",
            "(max|Pn[D*]| = ", signif(diag$ratio, 3), " standard errors). ",
            "Estimates are returned as NA.", call. = FALSE)
  }

  final <- NULL
  if(identical(status, "converged") && isTRUE(spec$supports_bayes)) {
    psi_f <- spec$psi_from_state(problem, state)
    Minv_f <- normalizing_matrix(problem$Lm_fn, psi_f, state$beta,
                                 problem$design_matrix, state$Q, problem$p)
    clever_f <- spec$make_clever(problem, steps[[1]], state, psi_f, Minv_f)
    clever_f$offset <- spec$fluctuation_offset(problem, steps[[1]], state)

    final <- list(
      psi = psi_f,
      Minv = Minv_f,
      K_Q = calculate_K(problem$Lm_fn, psi_f, state$beta,
                        problem$design_matrix, Minv_f, problem$batched_beta),
      clever = clever_f,
      epsilon = eps_list[[1]] # For MCM initialization
    )
  }

  list(state = state, status = status,
       converged = identical(status, "converged"),
       iter = iter, steps = steps, diag = diag,
       solver_log = do.call(rbind, lapply(solver_log, as.data.frame)),
       final = final)
}

finalize_tmle <- function(problem, spec, fit) {
  n <- problem$n
  p <- problem$p

  if(!identical(fit$status, "converged")) {
    return(list(
      est = rep(NA_real_, p),
      se = rep(NA_real_, p),
      lower = rep(NA_real_, p),
      upper = rep(NA_real_, p),
      eif = matrix(NA_real_, nrow = n, ncol = p),
      converged = FALSE,
      status = fit$status,
      iter = fit$iter,
      solved = if(is.null(fit$diag)) NA_real_ else fit$diag$ratio,
      solver_log = fit$solver_log,
      design_invariant = problem$design_invariant
    ))
  }

  d <- fit$diag
  est <- as.numeric(d$est)
  se <- as.numeric(d$se)

  list(
    est = est,
    se = se,
    lower = est + stats::qnorm(0.025) * se,
    upper = est + stats::qnorm(0.975) * se,
    eif = d$eif,
    converged = TRUE,
    status = "converged",
    iter = fit$iter,
    solved = d$ratio,
    solver_log = fit$solver_log,
    design_invariant = problem$design_invariant
  )
}
