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
