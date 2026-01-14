bound <- function(x, min, max, epsilon) {
  x[x <= min + epsilon] <- epsilon
  x[x >= max - epsilon] <- epsilon
  x
}
