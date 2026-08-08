bound <- function(x, min, max, epsilon) {
  x[x <= min + epsilon] <- epsilon
  x[x >= max - epsilon] <- epsilon
  x
}

bound01 <- function(x, epsilon) {
  bound(x, 0, 1, epsilon)
}
