`%||%` <- function(x, y) if(is.null(x)) y else x

bound <- function(x, min, max, epsilon) {
  x[x <= min + epsilon] <- epsilon
  x[x >= max - epsilon] <- epsilon
  x
}

bound01 <- function(x, epsilon) {
  bound(x, 0, 1, epsilon)
}

#' Coerce to a float32 torch tensor, preserving shape
#' @noRd
as_float_tensor <- function(x) {
  if(inherits(x, "torch_tensor")) return(x$to(dtype = torch::torch_float()))
  torch::torch_tensor(x, dtype = torch::torch_float())
}

#' Resolve fluctuation = "auto" once outcome_type is known
#' @noRd
resolve_fluctuation <- function(tmle, outcome_type) {
  if(identical(tmle$fluctuation, "auto")) {
    tmle$fluctuation <- if(identical(outcome_type, "binomial")) "logistic" else "linear"
  }
  if(identical(tmle$fluctuation, "logistic") && identical(outcome_type, "continuous")) {
    warning("A logistic fluctuation with `outcome_type` = \"continuous\"` ",
            "requires Y in [0, 1]; rescale Y first or use `fluctuation = \"linear\"`.", call. = FALSE)
  }
  tmle$linear <- identical(tmle$fluctuation, "linear")
  tmle
}
