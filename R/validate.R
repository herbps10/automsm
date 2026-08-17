# ----- Custom assertions -----

assert_sl_library <- function(x, name) {
  if(is.character(x)) {
    checkmate::assert_character(x, min.len = 1L, any.missing = FALSE, .var.name = name)
  }
  else if(is.list(x)) {
    checkmate::assert_list(x, min.len = 1L, .var.name = name)
    for(el in x) {
      checkmate::assert_character(el, min.len = 1L, any.missing = FALSE, .var.name = paste0(name, " element"))
    }
  }
  else {
    stop("`", name, "` must be a character vector of SuperLearner wrapper names, or a list of such vectors.", call. = FALSE)
  }
  invisible(TRUE)
}

# ----- Validations -----

#' Validate arguments for longitudinal  NP-MSM estimators
#' @noRd
validate_longitudinal_arguments <- function(
    data, Ls, As, Y, formula, summary_measures,
    outcome_type, loss, working_model, p,
    nuisance_estimates) {
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)
  checkmate::assert_character(As, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(As, choices = names(data))
  checkmate::assert_list(Ls, len = length(As))
  for(Lt in Ls) checkmate::assert_subset(Lt, choices = names(data))
  checkmate::assert_string(Y)
  checkmate::assert_choice(Y, choices = names(data))
  checkmate::assert_formula(formula)
  checkmate::assert_function(summary_measures, null.ok = TRUE)
  checkmate::assert_choice(outcome_type, choices = c("continuous", "binomial"))
  checkmate::assert_function(loss)
  checkmate::assert_function(working_model)

  checkmate::assert(
    checkmate::check_null(nuisance_estimates),
    checkmate::check_list(nuisance_estimates),
    combine = "or"
  )
}

#' Validate common arguments for single-time-point NP-MSM estimators
#'
#' @noRd
validate_msm_arguments <- function(
    data, X, A, Y, formula,
    outcome_type, loss, working_model, p,
    nuisance_estimates) {
  # Argument checks
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 1)
  checkmate::assert_string(A)
  checkmate::assert_subset(A, choices = names(data))
  checkmate::assert_character(X, min.len = 1, any.missing = TRUE, unique = TRUE)
  checkmate::assert_subset(X, choices = names(data))
  checkmate::assert_string(Y)
  checkmate::assert_choice(Y, choices = names(data))
  checkmate::assert_formula(formula)
  checkmate::assert_choice(outcome_type, choices = c("continuous", "binomial"))
  checkmate::assert_function(loss)
  checkmate::assert_function(working_model)
  checkmate::assert_count(p, positive = TRUE, null.ok = TRUE)

  checkmate::assert(
    checkmate::check_null(nuisance_estimates),
    checkmate::check_list(nuisance_estimates),
    combine = "or"
  )

  TRUE
}
