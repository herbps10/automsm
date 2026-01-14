#' @importFrom torch torch_tensor optim_lbfgs
#' @noRd
tmle_mle <- function(p, fluctuation_model, ...) {
  epsilon <- torch::torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(epsilon)
  for(iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      target <- fluctuation_model(epsilon, ...)
      #cat(glue::glue("TMLE iteration {iter} target: {as.numeric(target)}\n\n"))
      target$backward(retain_graph = TRUE)
      target
    })
  }
  optimizer$zero_grad()
  epsilon
}
