#' Raw clever-covariate directions
#'
#' Weighted parameterization
#' No normalizing matrix, and no inverse-propensity weight.
#' In this parameterization, W / prod(pi) goes in the GLM weights.
#' @noRd
nabla_Ldot <- function(problem, psi, beta) {
  blocks <- batched_NablaLdot(problem$Lm_fn, psi, beta,
                              problem$design_matrix, problem$p, problem$K)
  torch::torch_stack(blocks, dim = 1)$permute(c(2, 3, 1))
}

#' @noRd
fluctuate_ice <- function(mu, lin, control) {
  if(control$linear) return(mu + lin)
  b <- control$clamp %||% 1e-4
  pmin(pmax(stats::plogis(stats::qlogis(mu) + lin), b), 1 - b)
}


#' @noRd
lin_slice <- function(Xa, idx, eps) {
  m <- length(idx)
  K <- dim(Xa)[2]
  matrix(matrix(Xa[idx, , , drop = FALSE], nrow = m * K) %*% eps, nrow = m, ncol = K)
}

#' Telescoping sequential residual from out-of-fold array
#' @noRd
ice_delta <- function(problem, nodes) {
  n <- problem$n
  K <- problem$K
  tau <- problem$tau
  resid <- array(0, dim = c(K, n, tau))
  for(t in seq_len(tau)) {
    resid[, , t] <- nodes[, , t + 1L] - nodes[, , t]
  }
  resid <- resid / problem$aux$pi_cumprod * problem$aux$W
  stopifnot(all(is.finite(resid)))
  t(apply(resid, c(1, 2), sum))
}

#' One backward pass
#' @noRd
ice_backward_pass <- function(problem, control, beta, psi) {
  engine <- problem$nuisance_engine
  folds <- engine$folds
  V <- length(folds)

  n <- problem$n
  K <- problem$K
  p <- problem$p
  tau <- problem$tau
  family <- fluctuation_family(control$linear)

  Xa <- as.array(nabla_Ldot(problem, psi, beta)$detach())
  stopifnot(identical(dim(Xa), c(n, K, p)))

  # column-major throughout: (i, j) -> i + n * (j - 1)
  Xflat <- matrix(Xa, nrow = n * K)
  W <- problem$aux$W
  pic <- problem$aux$pi_cumprod
  Yv <- as.numeric(problem$Yt)

  st_train <- lapply(folds, function(f) {
    matrix(Yv[f$training_set], nrow = length(f$training_set), ncol = K)
  })

  st_valid <- lapply(folds, function(f) {
    matrix(Yv[f$validation_set], nrow = length(f$validation_set), ncol = K)
  })

  nodes <- array(NA_real_, dim = c(K, n, tau + 1L))
  for(j in seq_len(K)) {
    nodes[j, , tau + 1L] <- Yv
  }
  node_log <- vector("list", tau)

  for(t in tau:1) {
    # Fold-local regressions against current pseudo-outcome
    nd <- lapply(seq_len(V), function(v) engine$refit_node(t, st_train[[v]], v))

    # Stitch together pooled out-of-fold offset and target
    offset <- target <- matrix(NA_real_, n, K)
    for(v in seq_len(V)) {
      valid <- folds[[v]]$validation_set
      offset[valid, ] <- nd[[v]]$valid
      target[valid, ] <- st_valid[[v]]
    }
    stopifnot(!anyNA(offset), !anyNA(target))

    # Fit pooled epsilon, from the out-of-fold predictions
    wt <- as.vector(t(W[, , t] / pic[, , t]))
    keep <- wt > 0
    glm_t <- list(
      X = Xflat[keep, , drop = FALSE],
      offset = if(control$linear) as.vector(offset)[keep] else stats::qlogis(as.vector(offset)[keep]),
      target = as.vector(target)[keep],
      weights = wt[keep]
    )

    solution <- solve_fluctuation_glm(p, glm_t, family, control, paste0("node t=", t))
    eps <- solution$epsilon

    # Apply epsilon to every fold's training and validation slices
    for(v in seq_len(V)) {
      train <- folds[[v]]$training_set
      valid <- folds[[v]]$validation_set

      st_train[[v]] <- fluctuate_ice(nd[[v]]$train, lin_slice(Xa, train, eps), control)
      st_valid[[v]] <- fluctuate_ice(nd[[v]]$valid, lin_slice(Xa, valid, eps), control)
    }

    mu_t <- matrix(NA_real_, n, K)
    for(v in seq_len(V)) {
      mu_t[folds[[v]]$validation_set, ] <- st_valid[[v]]
    }
    nodes[, , t] <- t(mu_t)

    node_log[[t]] <- c(list(t = t, nnz = sum(keep)), solution$log)
  }

  list(nodes = nodes, psi = t(nodes[, , 1L]), node_log = do.call(rbind, lapply(node_log, as.data.frame)))
}

#' @noRd
eif_diagnostics_ice <- function(problem, nodes, psi_star, beta) {
  psi <- as_float_tensor(psi_star)$detach()$clone()$requires_grad_(TRUE)
  Delta <- as_float_tensor(ice_delta(problem, nodes))
  eifm <- eif(problem$Lm_fn, psi, beta, problem$design_matrix, problem$Q0, Delta, problem$p, problem$batched_beta)
  se <- apply(eifm, 2, stats::sd) / sqrt(problem$n)
  list(
    est = beta,
    eif = eifm,
    se = se,
    ratio = max(abs(colMeans(eifm)) / pmax(se, .Machine$double.eps))
  )
}

#' @noRd
ice_beta <- function(problem, psi_star) {
  psi <- as_float_tensor(psi_star)$detach()$clone()$requires_grad_(TRUE)
  b <- B(problem$Lm_fn, psi, problem$design_matrix, problem$Q0, problem$p)$detach()$clone()
  b$requires_grad_(TRUE)
  b
}

run_ice_tmle <- function(problem, spec, control, state0 = NULL) {
  engine <- problem$nuisance_engine

  if(!isTRUE(engine$supports_refit)) {
    stop("The longitudinal TMLE re-estimates the sequential regressions against ",
         "the targeted pseudo-outcomes, so it requires a nuisance engine that ",
         "supports refitting. Supply `nuisance = nuisance_control(...)` rather ",
         "than `nuisance_estimates`.", call. = FALSE)
  }

  if(identical(control$criterion, "epsilon")) {
    stop("`criterion = \"epsilon\"` is not available for the longitudinal estimand. Use \"eif\" (default) or \"beta\".", call. = FALSE)
  }


  maxit <- if(isTRUE(problem$design_invariant)) 1L else control$maxiter

  psi_cur <- t(problem$nuisance_estimates$mu[, , 1L])
  beta <- state0$beta
  beta_prev <- as.numeric(beta)
  logs <- list()

  for(k in seq_len(maxit)) {
    out <- ice_backward_pass(problem, control, beta,
                             as_float_tensor(psi_cur)$detach()$clone()$requires_grad_(TRUE))
    beta <- ice_beta(problem, out$psi)
    diag <- eif_diagnostics_ice(problem, out$nodes, out$psi, beta)

    logs[[k]] <- cbind(pass = k, out = out$node_log)

    dbeta <- max(abs(as.numeric(beta) - beta_prev))
    beta_prev <- as.numeric(beta)
    psi_cur <- out$psi
    done <- switch(control$criterion, eif = diag$ratio < control$eif_tol, beta = dbeta < (control$tol %||% 1e-6))
    if(done) break
  }

  status <- if(all(out$node_log$improved) && diag$ratio < control$eif_tol) {
    "converged"
  }
  else {
    "unsolved"
  }

  if(!identical(status, "converged")) {
    warning("Longitudinal TMLE did not solve the estimating equation ",
            "(max|Pn[D*]| = ", signif(diag$ratio, 3), " standard errors). ",
            "Estimates are returned as NA.", call. = FALSE)
  }

  list(
    state = list(nodes = out$nodes, beta = beta, Q = problem$Q0),
    status = status, converged = identical(status, "converged"),
    design_invariant = problem$design_invariant,
    iter = 1L, diag = diag, solver_log = out$node_log, final = NULL
  )

}
