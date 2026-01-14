#' Squared-error loss function constructor.
#' @importFrom torch nn_mse_loss
#' @param input input vector
#' @param target target vector
#' @details This is a wrapper of the torch::nn_mse_loss function (with reduction = "none").
#' @seealso [torch::nn_mse_loss()]
#' @export
loss_squared_error <- torch::nn_mse_loss(reduction = "none")

#' Smooth L1 loss function constructor.
#' @importFrom torch nn_smooth_l1_loss
#' @param input input vector
#' @param target target vector
#' @seealso [torch::nn_smooth_l1_loss()]
#' @details This is a wrapper of the torch::nn_smooth_l1_loss function (with reduction = "none").
#' @export
loss_smooth_l1 <- torch::nn_smooth_l1_loss(reduction = "none")
