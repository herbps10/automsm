autograd_gradient <- function(f, x) {
  e <- torch::torch_tensor(x, requires_grad = TRUE)
  as.numeric(torch::autograd_grad(f(e), e)[[1]])
}

hessian_from_gradient <- function(f, x, h = 1e-3) {
  p <- length(x)
  H <- matrix(0, p, p)
  for(j in seq_len(p)) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    H[, j] <- (autograd_gradient(f, xp) - autograd_gradient(f, xm)) / (2 * h)
  }
  attr(H, "asymmetry") <- max(abs(H - t(H))) / max(abs(H))
  (H + t(H)) / 2
}

Lhat_theory <- function(it) {
  st <- it$fit$state
  fn <- it$fit$final
  clever <- as.matrix(fn$clever$obs)
  K <- as.matrix(fn$K_Q)
  mu <- as.numeric(st$mu$obs)
  Q <- as.numeric(st$Q)
  n <- nrow(clever)

  A <- crossprod(as.matrix(clever) * sqrt(mu * (1 - mu)))

  Kbar <- colSums(Q * K)

  B <- n * (crossprod(K * sqrt(Q)) - tcrossprod(Kbar))
  A + B
}

bayes_loglik_at <- function(problem, spec, fit, condvar) {
  fn <- fit$final
  function(e) {
    st <- spec$apply_update(problem, list(id = 1L, fluctuate_Q = TRUE),
                            fit$state, e, fn$clever, fn$K_Q)
    clever_b <- scale_bayes_clever(problem, spec, st, fn$clever)
    spec$bayes_loglik(e, problem, fit$state, clever_b, st$Q, condvar)
  }
}
