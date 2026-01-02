

tmle_mle <- \(p, fluctuation_model, ...) {
  epsilon <- torch::torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(epsilon)
  for(iter in 1:5) {
    optimizer$step(\() {
      optimizer$zero_grad()
      target <- fluctuation_model(epsilon, ...)
      target$backward(retain_graph = TRUE)
      target
    })
  }
  optimizer$zero_grad()
  epsilon
}
