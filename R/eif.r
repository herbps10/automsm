#' Convenience notation for writing composition of loss function and working model
#' @noRd
Lm <- function(loss, working_model) {
  function(t, beta, X) {
    pred <- working_model(beta, X)
    l <- loss(working_model(beta, X), t)
    if(length(dim(l)) > 1) l$sum(2) else l
  }
}

#' Gradient of beta -> Lm(t, beta, X)
#' @importFrom torch autograd_grad
#' @noRd
dL <- function(Lm, t, beta, X, weight = 1) {
  beta$grad <- NULL
  l <- Lm(t, beta, X)$mul(weight)$sum()
  torch::autograd_grad(l, list(beta), create_graph = TRUE, retain_graph = TRUE)
}

#' @importFrom purrr map
#' @importFrom torch autograd_grad torch_cat
#' @noRd
dL_dbeta <- function(Lm, t, beta, X, p) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- purrr::map(1:length(l), function(index) {
    torch::autograd_grad(
      l[index],
      list(beta),
      retain_graph = TRUE,
      create_graph = TRUE
    )[[1]]$reshape(c(p, 1))
  })
  torch::torch_cat(r, dim = 2)
}

#' @importFrom torch torch_reshape autograd_grad torch_cat
#' @importFrom purrr map
#' @noRd
ddL <- function(Lm, t, beta, X, p, weight = 1) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X, weight = 1)[[1]]
  r <- map(1:length(l), function(index) {
    torch::torch_reshape(
      torch::autograd_grad(
        l[index],
        list(beta),
        retain_graph = TRUE,
        create_graph = TRUE
      )[[1]],
      shape = c(p, 1)
    )
  })
  torch::torch_cat(r, 2)
}

#' @importFrom purrr map
#' @importFrom torch torch_reshape autograd_grad torch_cat
#' @noRd
grad_dL <- function(Lm, t, beta, X) {
  k <- dim(X)[2]
  t$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- map(1:length(l), function(index) {
    torch::torch_reshape(
      torch::autograd_grad(
        l[index],
        list(t),
        retain_graph = TRUE,
        create_graph = TRUE,
        allow_unused = TRUE
      )[[1]],
      shape = c(1, k)
    )
  })
  torch::torch_cat(r)
}

#' @importFrom torch optim_lbfgs torch_tensor
#' @noRd
B <- function(Lm, psi, design_matrix, Q, p = NULL) {
  if(is.null(p)) p <- rev(dim(design_matrix))[1]
  beta <- torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(beta)
  for (iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      value <- Lm(psi, beta, design_matrix)$mul(Q)$sum()
      #cat(glue::glue("Iteration: {iter} value: {as.numeric(value)} \n\n"))
      value$backward(retain_graph = TRUE)
      value
    })
    optimizer$zero_grad()
  }
  beta
}

#' @noRd
objective <- function(Lm, psi, Q, design_matrix, beta) {
  dL(Lm, psi, beta, design_matrix, weight = Q)[[1]]
}

#' @importFrom torch torch_cat autograd_grad
#' @importFrom purrr map
#' @noRd
dobjective_dpsi <- function(Lm, psi, Q, design_matrix, beta, p = NULL) {
  n <- dim(design_matrix)[1]
  k <- dim(design_matrix)[2]
  if(is.null(p)) p <- dim(design_matrix)[3]

  obj <- objective(Lm, psi, Q, design_matrix, beta)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], psi, retain_graph = TRUE)[[1]]$reshape(c(
      1,
      n,
      k
    ))
  })
  return(torch::torch_cat(r)$transpose(1, 3))
}

#' @importFrom torch autograd_grad torch_cat
#' @importFrom purrr map
#' @noRd
dobjective_dQ <- function(Lm, psi, Q, design_matrix, beta, p = NULL) {
  if(is.null(p)) p <- rev(dim(design_matrix))[1]
  n <- dim(design_matrix)[1]

  obj <- objective(Lm, psi, Q, design_matrix, beta)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], Q, retain_graph = TRUE)[[1]]$reshape(c(n, 1))
  })
  torch::torch_cat(r, dim = 2)
}

#' @importFrom purrr map
#' @importFrom torch autograd_grad torch_cat
#' @noRd
dobjective_dbeta <- function(Lm, psi, Q, design_matrix, beta, p =  NULL) {
  if(is.null(p)) p <- dim(design_matrix)[3]
  obj <- objective(Lm, psi, Q, design_matrix, beta)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], beta, retain_graph = TRUE)[[1]]$reshape(c(p, 1))
  })
  torch::torch_cat(r, dim = 2)
}

dB_dpsi <- function(Lm, psi, Q, design_matrix, beta) {
  -dobjective_dpsi(
    Lm,
    psi,
    Q,
    design_matrix,
    beta
  )$matmul(torch::linalg_inv(dobjective_dbeta(Lm, psi, Q, design_matrix, beta)))
}

dB_dQ <- function(Lm, psi, Q, design_matrix, beta) {
  -dobjective_dQ(
    Lm,
    psi,
    Q,
    design_matrix,
    beta
  )$matmul(torch::linalg_inv(dobjective_dbeta(Lm, psi, Q, design_matrix, beta)))
}

#' @importFrom torch torch_tensor
#' @noRd
normalizing_matrix <- function(Lm, psi, beta, design_matrix, Q, p = NULL) {
  if(is.null(p)) p <- dim(design_matrix)[3]
  n <- dim(design_matrix)[1]
  M <- matrix(0, p, p)
  M <- ddL(Lm, psi, beta, design_matrix, p, Q) / n
  Minv <- solve(M)
  torch::torch_tensor(Minv)
}

#' Per-observation gradient of the loss wrt beta (the \dot{L}_m / D_2 term)
#' Returns an n x p tensor whose i-th row is dL_m/dbeta for observation i
#' @noRd
batched_dL_dbeta <- function(Lm, psi, beta, design_matrix, p) {
  n <- dim(design_matrix)[1]
  Lm_vec <- Lm(psi, beta, design_matrix)
  D2 <- torch::torch_cat(
    purrr::map(1:n, function(i) {
      beta$grad <- NULL
      sel <- torch::torch_zeros(n)
      sel[i] <- 1
      torch::autograd_grad(
        outputs = Lm_vec,
        inputs = beta,
        grad_outputs = sel,
        retain_graph = TRUE,
        create_graph = FALSE
      )[[1]]$reshape(c(1, p))
    }),
    dim = 1
  )
}

#' Per-observation NablaLdot blocks (the grad-of-loss-gradient-wrt-psi term)
#' Returns a list of p tensors, each n x k, with NULL gradients zero-filled
#' @noRd
batched_NablaLdot <- function(Lm, psi, beta, design_matrix, p, k) {
  n <- dim(design_matrix)[1]
  obj <- dL(Lm, psi, beta, design_matrix)[[1]]
  purrr::map(1:p, function(j) {
    g <- torch::autograd_grad(
      obj[j],
      psi,
      retain_graph = TRUE,
      allow_unused = TRUE
    )[[1]]

    if(is.null(g)) {
      g <- torch::torch_zeros(c(n, k))
    }
    else {
      g <- g$reshape(c(n, k))
    }
    g
  })
}

#' @noRd
eif <- function(Lm, psi, beta, design_matrix, Q, Delta, p = NULL) {
  n <- dim(design_matrix)[1]
  k <- dim(design_matrix)[2]
  if(is.null(p)) p <- dim(design_matrix)[3]

  Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q, p)

  D2 <- batched_dL_dbeta(Lm, psi, beta, design_matrix, p)

  # D1 = NablaLdot * Delta
  grad_blocks <- batched_NablaLdot(Lm, psi, beta, design_matrix, p, k)

  D1 <- -torch::torch_cat(
    purrr::map(grad_blocks, function(gb) {
      (gb * Delta)$sum(dim = 2)$reshape(c(n, 1))
    }),
    dim = 2
  )

  eif <- as.matrix((D1 + D2)$matmul(-Minv))
  return(eif)
}
