onestep <- function(Lm, psi, beta, design_matrix, Q, Delta) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]

  onestep_eif   <- eif(Lm, psi, beta, design_matrix, Q, Delta)
  onestep_est   <- beta + colMeans(onestep_eif)
  onestep_se    <- apply(onestep_eif, 2, sd) / sqrt(n)
  onestep_lower <- onestep_est + qnorm(0.025) * onestep_se
  onestep_upper <- onestep_est + qnorm(0.975) * onestep_se
  list(
    est = onestep_est,
    se = onestep_se,
    eif = onestep_eif,
    lower = onestep_lower,
    upper = onestep_upper
  )
}
