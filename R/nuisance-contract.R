#' @noRd
nuisance_field <- function(name, dim, required = TRUE, lower = NULL, upper = NULL, severity = "error") {
  list(name = name, dim = as.integer(dim), required = required, lower = lower, upper = upper, severity = severity)
}

#' @noRd
nuisance_problem <- function(severity, ...) {
  if(identical(severity, "warning")) {
    warning(..., call. = FALSE)
  }
  else {
    stop(..., call. = FALSE)
  }
  invisible(NULL)
}

#' Validate a user-supplied or fitted nuisance list against the estimand contract
#' @noRd
validate_nuisance <- function(problem, spec, bayes_enabled = FALSE) {
  nu <- problem$nuisance
  checkmate::assert_list(nu, names = "named", .var.name = "nuisance")
  contract <- spec$nuisance_contract(problem, bayes_enabled)

  for(f in contract$fields) {
    x <- nu[[f$name]]
    if(is.null(x)) {
      if(f$required) {
        stop("nuisance$", f$name, " is required for estimand '", problem$estimand,
             "'", if(bayes_enabled) " with bayes = TRUE" else "", ".", call. = FALSE)
      }
      next
    }

    got <- if(is.null(dim(x))) length(x) else dim(x)
    if(!identical(as.integer(got), f$dim)) {
      stop("nuisance$", f$name, " has shape(", paste(got, collapse = ", "),
           "); expected (", paste(f$dim, collapse = ", "), ").", call. = FALSE)
    }

    v <- as.numeric(x)
    if(anyNA(v) || !all(is.finite(v))) {
      stop("nuisance$", f$name, " contains missing or non-finite values.", call. = FALSE)
    }
    if(!is.null(f$lower) && min(v) < f$lower) {
      nuisance_problem(f$severity, "nuisance$", f$name, " must be > ", f$lower, "; observed minimum ", signif(min(v), 4), ".")
    }
    if(!is.null(f$upper) && max(v) > f$upper) {
      nuisance_problem(f$severity, "nuisance$", f$name, " must be < ", f$upper, "; observed maximum ", signif(max(v), 4), ".")
    }
  }

  for(check in contract$checks) {
    msg <- check(nu, problem)
    if(!isTRUE(msg)) nuisance_problem(attr(msg, "severity") %||% "error", msg)
  }

  invisible(TRUE)
}

#' @noRd
outcome_bounds <- function(problem) {
  if(identical(problem$outcome_type, "binomial")) list(lo = 0, hi = 1)
  else list(lo = NULL, hi = NULL)
}
