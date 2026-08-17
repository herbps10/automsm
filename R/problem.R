build_design_tensor <- function(formula, data, K, mutate = NULL, p = NULL) {
  n <- nrow(data)
  if(is.null(mutate)) {
    mat <- stats::model.matrix(formula, data = data)
    d <- ncol(mat)
    dm <- torch::torch_tensor(mat)$reshape(c(n, 1L, d))
    if(K > 1L) dm <- dm$expand(c(n, K, d))
    return(list(
      design_matrix = dm,
      terms = colnames(mat),
      d = d,
      p = p %||% d
    ))
  }
  mat1 <- stats::model.matrix(formula, data = mutate(data, 1L))
  d <- ncol(mat1)
  dm <- torch::torch_zeros(c(n, K, d))
  dm[, 1, ] <- mat1
  if(K > 1L) {
    for(k in 2:K) {
      dm[, k, ] <- stats::model.matrix(formula, data = mutate(data, k))
    }
  }

  list(
    design_matrix = dm,
    terms = colnames(mat1),
    d = d,
    p = p %||% d
  )
}

validate_working_model_dim <- function(working_model, design_matrix, p) {
  out <- try(working_model(torch::torch_zeros(p), design_matrix), silent = TRUE)

  if(inherits(out, "try-error")) {
    stop("working_model() failed with a length-", p, " beta. Supply `p` ",
         "matching the working model's coefficient count.", call. = FALSE)
  }

  if(!identical(dim(out), dim(design_matrix)[1:2])) {
    stop("working_model() returned shape(", paste(dim(out), collapse = ", "),
         "); expected (", paste0(dim(design_matrix[1:2]), collapse = ", "), ").", call. = FALSE)
  }
  TRUE
}

#' @noRd
validate_design_conditioning <- function(design_matrix, terms, d, error_kappa = 1e12, warn_kappa = 1e7) {
  n <- dim(design_matrix)[1]
  K <- dim(design_matrix)[2]
  X <- as.matrix(design_matrix$reshape(c(n * K, d)))

  zero <- which(apply(X, 2, function(cc) max(abs(cc))) < 1e-10)

  if(length(zero)) {
    stop("Working-model design column(s) ", paste(zero, collapse = ", "),
         " (", paste(terms[zero], collapse = ", "), ") are identically zero. ",
    call. = FALSE)
  }

  k <- kappa(crossprod(X), exact = TRUE)
  if(!is.finite(k) || k > error_kappa) {
    stop("Working-model design is numerically rank-deficient ",
         "(condition number ", signif(k, 3), ", rank ", qr(X)$rank, " of ", d,
         "). The normalizing matrix M is then singular, so B(P) is not ",
         "uniquely defined and beta is not identified. Drop collinear terms or respecify `formula`.")
  }

  if(k > warn_kappa) {
    warning("Working-model design is ill-conditioned (condition number ", signif(k, 3), "). Estimates and standard errors may be unstable.", call. = FALSE)
  }

  invisible(k)
}

probe_batched_beta <- function(working_model, design_matrix, p) {
  m <- min(dim(design_matrix[1]), 8L)
  Xs <- design_matrix[1:m, , , drop = FALSE]
  b <- torch::torch_randn(p)
  one <- try(working_model(b, Xs), silent = TRUE)
  if(inherits(one, "try-error")) return(FALSE)
  Bb <- b$unsqueeze(1)$expand(c(m, p))$clone()
  many <- try(working_model(Bb, Xs), silent = TRUE)
  if(inherits(many, "try-error") || !identical(dim(many), dim(one))) return(FALSE)
  isTRUE(as.numeric((many - one)$abs()$max()) < 1e-5)
}

new_msm_problem <- function(estimand, K, d, p, tau = 1L,
                            design_matrix, Q0, Yt, Lm_fn,
                            loss, working_model, formula, terms,
                            outcome_type, nuisance_estimates, aux = list()) {
  n <- dim(design_matrix)[1]
  stopifnot(identical(dim(design_matrix), c(n, as.integer(K), as.integer(d))))
  validate_working_model_dim(working_model, design_matrix, p)
  design_kappa <- validate_design_conditioning(design_matrix, terms, d)
  batched_beta <- probe_batched_beta(working_model, design_matrix, p)
  structure(
    list(estimand = estimand, n = n, K = as.integer(K), d = as.integer(d),
         p = as.integer(p), tau = as.integer(tau),
         design_matrix = design_matrix, Q0 = Q0, Yt = as_float_tensor(Yt),
         Lm_fn = Lm_fn, loss = loss, working_model = working_model,
         formula = formula, terms = terms, outcome_type = outcome_type,
         batched_beta = batched_beta,
         design_kappa = design_kappa,
         nuisance_estimates = nuisance_estimates, aux = aux),
    class = "msm_problem"
  )
}
