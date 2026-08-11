#' Squared-error loss
#'
#' Elementwise squared-error loss, \eqn{L(a, b) = (a - b)^2}, suitable
#' for continuous or bounded outcomes. This is a thin wrapper around
#' \code{torch::nn_mse_loss} with \code{reduction = "none"}, so
#' it returns per-element loss rather than a reduced scalar.
#'
#' @details The returned object is a \code{torch} loss module. It is called
#' as \code{loss_squared_error(prediction, target)}, where both arguments are
#' \code{torch} tensors of matching shape. \code{nn_mse_loss} is symmetric in
#' its two arguments.
#'
#' @usage loss_squared_error
#'
#' @importFrom torch nn_mse_loss
#' @seealso [torch::nn_mse_loss()], [loss_weighted_sum()]
#' @export
loss_squared_error <- torch::nn_mse_loss(reduction = "none")

#' Binary cross-entropy loss
#'
#' Elementwise binary cross-entropy loss suitable for binary outcomes.
#' This is a thin wrapper around \code{torch::nn_bce_loss}
#' with \code{reduction = "none"}, so it returns per-element loss rather
#' than a reduced scalar.
#'
#' @details The returned object is a \code{torch} loss module. It is called
#' as \code{loss_cross_entropy(prediction, target)}, where
#' \code{prediction} is a tensor in \eqn{[0, 1]} and \code{target} is a tensor of
#' outcomes of matching shape. The argument order is meaningful:
#' the first argument must be the predictions.
#'
#' @usage loss_cross_entropy
#'
#' @seealso [torch::nn_bce_loss()], [loss_weighted_sum()]
#' @importFrom torch nn_bce_loss
#' @export
loss_cross_entropy <- torch::nn_bce_loss(reduction = "none")

#' Binary cross-entropy loss (logit parameterization)
#'
#' Elementwise binary cross-entropy loss where the prediction argument is on
#' the logit scale, suitable for binary outcomes. This is a thin wrapper around
#' \code{torch::nn_bce_with_logits_loss} with \code{reduction = "none"}, so
#' it returns per-element loss rather than a reduced scalar.
#'
#' @details The returned object is a \code{torch} loss module. It is called
#' as \code{loss_cross_entropy_logit(prediction, target)}, where
#' \code{prediction} is a tensor of logits and \code{target} is a tensor of
#' outcomes of matching shape. The argument order is meaningful:
#' the first argument must be the logits.
#'
#' @usage loss_cross_entropy_logit
#'
#' @seealso [torch::nn_bce_with_logits_loss()], [loss_weighted_sum()]
#' @importFrom torch nn_bce_with_logits_loss
#' @export
loss_cross_entropy_logit <- torch::nn_bce_with_logits_loss(reduction = "none")

#' Smooth L1 loss
#'
#' Elementwise smooth L1 (Huber-type) loss, a smooth approximation to the
#' absolute-error loss that behaves quadratically for small residuals and
#' linearly for large ones. This is a thin wrapper aoround
#' \code{torch::nn_smooth_l1_loss} with \code{reduction = "none"}, so it
#' returns per-element losses rather than a reduced scalar.
#'
#' @details The returned object is a \code{torch} loss module. It is called as
#' \code{loss_smooth_l1(prediction, target)}, where both arguments are
#' \code{torch} tensors of matching shape.
#'
#' @usage loss_smooth_l1
#'
#' @seealso [torch::nn_smooth_l1_loss()], [loss_weighted_sum()]
#' @importFrom torch nn_smooth_l1_loss
#' @export
loss_smooth_l1 <- torch::nn_smooth_l1_loss(reduction = "none")

#' Construct a weighted-sum loss over target-functional components
#'
#' Builds a multi-dimensional loss function by summing a one-dimensional
#' (elementwise) loss over each of the \eqn{K} components of a multi-dimensional
#' target functional, optionally weighting each component. This is the
#' natural way to form a loss for target functionals with \eqn{K > 1}, such as
#' the conditional mean under a longitudinal treatment, and lets
#' estimation functions such as \code{longitudinal_dose_response} remain
#' agnostic to how the multi-dimensional loss is assembled.
#'
#' @param loss A one-dimensional (elementwise) loss, such as [loss_squared_error()]
#'   or [loss_cross_entropy_logit()], called as \code{loss(prediction, target)}.
#' @param weights A weight function mapping a component index \eqn{j \in \{1, \dots, K\}}
#'   to a non-negative weight. The default (\code{function(x) 1})
#'   assigns equal weight to every element. To weight by empirical support,
#'   for example, supply a function that returns the proportion of observations
#'   following component \eqn{j}.
#'
#'   Note that, unless the working model is assumed correctly specified,
#'   the choice of weight function changes the target parameter being
#'   estimated (see Petersen et al., 2014). Also note that, when `weights` is
#'   estimated from the data, the influence-curve based inference is valid
#'   for the estimand defined by the estimated weights, and does not account
#'   for uncertainty in the estimation of `weights` itself.
#' @returns A loss function with signature \code{function(prediction, target)}
#'   returning a length-\eqn{n} \code{torch} tensor of per-observation summed
#'   losses, where \code{prediction} and \code{target} are each \eqn{n \times K}
#'   tensors (the assembled working-model outputs and target-functional values).
#' @seealso [loss_squared_error()], [loss_cross_entropy_logit()], [loss_smooth_l1()]
#' @importFrom torch torch_tensor
#' @export
loss_weighted_sum <- function(loss, weights = function(x) 1) {
  function(prediction, target) {
    n <- dim(prediction)[1]
    K <- dim(prediction)[2]
    total <- torch::torch_tensor(rep(0, n), requires_grad = TRUE)
    for(j in 1:K) {
      total <- total$add(loss(prediction[, j], target[, j])$mul(weights(j)))
    }
    total
  }
}
