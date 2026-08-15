test_that("calculate_K returns correct fluctuation covariate", {
  loss <- loss_squared_error
  working_model <- working_model_linear

  Minv <- torch::torch_tensor(matrix(c(1, 0.5, 0.5, 1), ncol = 2))$inverse()

  psi <- torch::torch_tensor(rep(1, 5e2))$reshape(c(500, 1))
  beta <- torch::torch_tensor(c(1, 1), requires_grad = TRUE)
  design_matrix <- torch::torch_tensor(array(1, dim = c(500, 1, 2)))

  K <- calculate_K(Lm(loss, working_model), psi, beta, design_matrix, Minv, TRUE)
  expect_equal(torch::as_array(K), matrix(1 + 1/3, ncol = 2, nrow = 5e2), tolerance = 1e-4)

  K <- calculate_K(Lm(loss, working_model), psi, beta, design_matrix, Minv, FALSE)
  expect_equal(torch::as_array(K), matrix(1 + 1/3, ncol = 2, nrow = 5e2), tolerance = 1e-4)
})
