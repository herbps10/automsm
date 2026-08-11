#' @importFrom stats qnorm
#' @noRd
onestep <- function(Lm, psi, beta, design_matrix, Q, Delta, p = NULL) {
  n <- dim(design_matrix)[1]
  if(is.null(p)) p <- dim(design_matrix)[3]

  onestep_eif <- eif(Lm, psi, beta, design_matrix, Q, Delta, p)
  onestep_est <- beta + colMeans(onestep_eif)
  onestep_se <- apply(onestep_eif, 2, sd) / sqrt(n)
  onestep_lower <- onestep_est + stats::qnorm(0.025) * onestep_se
  onestep_upper <- onestep_est + stats::qnorm(0.975) * onestep_se
  list(
    est = onestep_est,
    se = onestep_se,
    eif = onestep_eif,
    lower = onestep_lower,
    upper = onestep_upper
  )
}
