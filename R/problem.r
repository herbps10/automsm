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
         "); expected (", paste0(dim(design_matrix[1:2], collapse = ", "), ")."), call. = FALSE)
  }
  TRUE
}

new_msm_problem <- function(estimand, K, d, p, tau = 1L,
                            design_matrix, Q0, Yt, Lm_fn,
                            loss, working_model, formula, terms,
                            outcome_type, nuisance, aux = list()) {
  n <- dim(design_matrix)[1]
  stopifnot(identical(dim(design_matrix), c(n, as.integer(K), as.integer(d))))
  validate_working_model_dim(working_model, design_matrix, p)
  structure(
    list(estimand = estimand, n = n, K = as.integer(K), d = as.integer(d),
         p = as.integer(p), tau = as.integer(tau),
         design_matrix = design_matrix, Q0 = Q0, Yt = as_float_tensor(Yt),
         Lm_fn = Lm_fn, loss = loss, working_model = working_model,
         formula = formula, terms = terms, outcome_type = outcome_type,
         nuisance = nuisance, aux = aux),
    class = "msm_problem"
  )
}
