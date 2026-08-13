test_that("dL gives correct result for squared error loss with linear working model", {
  loss <- loss_squared_error
  working_model <- working_model_linear

  psi <- torch::torch_tensor(c(1, 1))
  beta <- torch::torch_tensor(c(1, 1), requires_grad = TRUE)
  design_matrix <- torch::torch_tensor(matrix(c(1, 1, 1, 1), ncol = 2))

  x <- dL(Lm(loss, working_model), psi, beta, design_matrix)[[1]]

  expect_equal(torch::as_array(x), c(4, 4))
})
