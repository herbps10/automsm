#' @export
loss_squared_error <- torch::nn_mse_loss(reduction = "none")

#' @export
loss_smooth_l1 <- torch::nn_smooth_l1_loss(reduction = "none")
