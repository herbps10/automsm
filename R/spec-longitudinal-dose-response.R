msm_spec_longitudinal_dose_response <- function(tmle_linear = TRUE) {
  new_msm_spec(
    estimand = "longitudinal_dose_response",
    supports_bayes = FALSE,
    driver = run_ice_tmle,

    init_state = function(problem) {
      list(
        nodes = problem$nuisance_estimates$mu,
        Q = problem$Q0,
        beta = NULL
      )
    },

    psi_from_state = function(problem, state) {
      psi <- as_float_tensor(t(state$nodes[, , 1L]))
      psi$detach()$clone()$requires_grad_(TRUE)
    },

    # Telescoping ICE residual:
    #  Delta[i, j] = sum_t W[j, i, t] (mu_{t+1} - mu_t)[j, i] / g_{0:t}[j, i]
    delta = function(problem, state) {
      n <- problem$n
      k <- problem$K
      tau <- problem$tau
      resid <- array(0, dim = c(k, n, tau))
      for(t in seq_len(tau)) {
        resid[, , t] <- state$nodes[, , t + 1L] - state$nodes[, , t]
      }
      resid <- resid / problem$aux$pi_cumprod * problem$aux$W
      stopifnot(all(is.finite(resid)))
      as_float_tensor(t(apply(resid, c(1, 2), sum))) # (n, k)
    },

    nuisance_contract = function(problem, bayes_enabled) {
      n <- problem$n
      K <- problem$K
      tau <- problem$tau
      b <- outcome_bounds(problem)
      list(
        fields = list(
          nuisance_field("pi", c(n, tau), lower = 0, upper = 1),
          nuisance_field("mu", c(K, n, tau + 1), lower = b$lo, upper = b$hi, severity = "warning")
        ),
        checks = list(
          function(nu, problem) {
            Y <- as.numeric(problem$Yt)
            ok <- vapply(seq_len(K), function(j) {
              max(abs(nu$mu[j, , problem$tau + 1L] - Y)) < 1e-6
            }, logical(1))
            if(all(ok)) TRUE else "nuisance$mu[, , tau + 1] must equal Y for every regime."
          }
        )
      )
    }
  )
}
