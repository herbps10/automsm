Lm <- function(loss, working_model) function(t, beta, X) loss(t, working_model(beta, X))

# Gradient of beta -> Lm(t, beta, X)
dL <- function(Lm, t, beta, X, weight = 1) {
  beta$grad <- NULL
  l <- Lm(t, beta, X)$mul(weight)$sum()
  torch::autograd_grad(l, list(beta), create_graph = TRUE, retain_graph = TRUE)
}

dL_dbeta <- function(Lm, t, beta, X, p) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- map(1:length(l), \(index)
      torch::autograd_grad(l[index], list(beta), retain_graph = TRUE, create_graph = TRUE)[[1]]$reshape(c(p, 1))
  )
  torch::torch_cat(r, dim = 2)
}

ddL <- function(Lm, t, beta, X, p, weight = 1) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X, weight = 1)[[1]]
  r <- map(1:length(l), \(index)
    torch::torch_reshape(torch::autograd_grad(l[index], list(beta), retain_graph = TRUE, create_graph = TRUE)[[1]], shape = c(p, 1))
  )
  torch::torch_cat(r, 2)
}

# ∇dL(t, beta, X) = ForwardDiff.jacobian(t -> dL(first(t), beta, X), [t])
grad_dL <- function(Lm, t, beta, X) {
  k <- 1
  if(length(dim(X)) == 3) k <- dim(X)[1]
  t$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- map(1:length(l), \(index) {
    torch::torch_reshape(torch::autograd_grad(l[index], list(t), retain_graph = TRUE, create_graph = TRUE, allow_unused = TRUE)[[1]], shape = c(1, k))
  })
  torch::torch_cat(r)
}

B <- function(Lm, psi, design_matrix, Q) {
  p <- rev(dim(design_matrix))[1]
  beta <- torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(beta)
  for(iter in 1:2) {
    optimizer$step(\() {
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

objective <- function(Lm, psi, Q, design_matrix, beta) {
  dL(Lm, psi, beta, design_matrix, weight = Q)[[1]]
}

dobjective_dpsi <- function(Lm, psi, Q, design_matrix, beta) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]

  if(length(dim(design_matrix)) == 3) {
    K <- dim(design_matrix)[1]
    obj <- objective(Lm, psi, Q, design_matrix, beta)
    r <- map(1:p, \(index) autograd_grad(obj[index], psi, retain_graph = TRUE)[[1]]$reshape(c(1, n, K)))
    return(torch::torch_cat(r)$transpose(1, 3))
  }
  else {
    obj <- objective(Lm, psi, Q, design_matrix, beta)
    r <- map(1:p, \(index) autograd_grad(obj[index], psi, retain_graph = TRUE)[[1]]$reshape(c(n, 1)))
    return(torch::torch_cat(r, dim = 2))
  }
}

dobjective_dQ <- function(Lm, psi, Q, design_matrix, beta) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]

  obj <- objective(Lm, psi, Q, design_matrix, beta)
  r <- map(1:p, \(index) autograd_grad(obj[index], Q, retain_graph = TRUE)[[1]]$reshape(c(n, 1)))
  torch::torch_cat(r, dim = 2)
}

dobjective_dbeta <- function(Lm, psi, Q, design_matrix, beta) {
  p <- rev(dim(design_matrix))[1]
  obj <- objective(Lm, psi, Q, design_matrix, beta)
  r <- map(1:p, \(index) autograd_grad(obj[index], beta, retain_graph = TRUE)[[1]]$reshape(c(p, 1)))
  torch::torch_cat(r, dim = 2)
}

dB_dpsi <- function(Lm, psi, Q, design_matrix, beta) {
  -dobjective_dpsi(Lm, psi, Q, design_matrix, beta)$matmul(torch::linalg_inv(dobjective_dbeta(Lm, psi, Q, design_matrix, beta)))
}

dB_dQ <- function(Lm, psi, Q, design_matrix, beta) {
  -dobjective_dQ(Lm, psi, Q, design_matrix, beta)$matmul(torch::linalg_inv(dobjective_dbeta(Lm, psi, Q, design_matrix, beta)))
}

normalizing_matrix <- function(Lm, psi, beta, design_matrix, Q) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]
  M <- matrix(0, p, p)
  #for(i in 1:n) {
  #  if(length(dim(design_matrix)) == 3) {
  #    M <- M + Q[i] * ddL(Lm, psi[i,, drop = FALSE], beta, design_matrix[, i, , drop=FALSE], p)
  #  }
  #  else {
  #    M <- M + Q[i] * ddL(Lm, psi[i], beta, design_matrix[i, ], p)
  #  }
  #}
  M <- ddL(Lm, psi, beta, design_matrix, p, Q) / n
  Minv <- solve(M)
  torch::torch_tensor(Minv)
}


eif <- function(Lm, psi, beta, design_matrix, Q, Delta) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]
  Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q)

  eif <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n) {
    if(length(dim(design_matrix)) == 3) {
      eif[i, ] <- -as.numeric(Minv$matmul(torch::torch_reshape(dL(Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE])[[1]], shape = c(p, 1)) + grad_dL(Lm, psi[i, drop = FALSE], beta, design_matrix[, i, drop = FALSE])$matmul(Delta[i,])$reshape(c(p, 1))))
    }
    else {
      eif[i, ] <- -as.numeric(Minv$matmul(torch::torch_reshape(dL(Lm, psi[i], beta, design_matrix[i, ])[[1]], shape = c(p, 1)) + grad_dL(Lm, psi[i], beta, design_matrix[i, ])) * Delta[i])
    }
  }
  eif
}
