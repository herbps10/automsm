Lm <- \(loss, working_model) \(t, beta, X) loss(t, working_model(beta, X))

# Gradient of beta -> Lm(t, beta, X)
dL <- \(Lm, t, beta, X) {
  beta$grad <- NULL
  l <- Lm(t, beta, X)
  torch::autograd_grad(l, list(beta), create_graph = TRUE, retain_graph = TRUE)
}

ddL <- \(Lm, t, beta, X, p) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- map(1:length(l), \(index)
    torch::torch_reshape(torch::autograd_grad(l[index], list(beta), retain_graph = TRUE, create_graph = TRUE)[[1]], shape = c(p, 1))
  )
  torch::torch_cat(r, 2)
}

# ∇dL(t, beta, X) = ForwardDiff.jacobian(t -> dL(first(t), beta, X), [t])
grad_dL <- \(Lm, t, beta, X) {
  k <- 1
  if(length(dim(X)) == 3) k <- dim(X)[1]
  t$grad <- NULL
  l <- dL(Lm, t, beta, X)[[1]]
  r <- map(1:length(l), \(index) {
    torch::torch_reshape(torch::autograd_grad(l[index], list(t), retain_graph = TRUE, create_graph = TRUE)[[1]], shape = c(1, k))
  })
  torch::torch_cat(r)
}


B <- \(Lm, psi, design_matrix, Q) {
  p <- rev(dim(design_matrix))[1]
  beta <- torch_tensor(rep(0, p), requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(beta)
  for(iter in 1:10) {
    optimizer$step(\() {
      optimizer$zero_grad()
      value <- Lm(psi, beta, design_matrix)$mul(Q)$sum()
      cat(glue::glue("Iteration: {iter} value: {as.numeric(value)} \n\n"))
      value$backward(retain_graph = TRUE)
      value
    })
  }
  optimizer$zero_grad()
  beta
}

normalizing_matrix <- \(Lm, psi, beta, design_matrix, Q) {
  p <- rev(dim(design_matrix))[1]
  n <- rev(dim(design_matrix))[2]
  M <- matrix(0, p, p)
  for(i in 1:n) {
    if(length(dim(design_matrix)) == 3) {
      M <- M + Q[i] * ddL(Lm, psi[i,, drop = FALSE], beta, design_matrix[, i, , drop=FALSE], p)
    }
    else {
      M <- M + Q[i] * ddL(Lm, psi[i], beta, design_matrix[i, ], p)
    }
  }
  Minv <- solve(M)
  torch::torch_tensor(Minv)
}


eif <- \(Lm, psi, beta, design_matrix, Q, Delta) {
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
