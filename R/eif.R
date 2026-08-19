# `p` is the working-model coefficient count and is always required.
# It is distinct from d = dim(design_matrix)[3], the design width.

#' Convenience notation for writing composition of loss function and working model
#' @noRd
Lm <- function(loss, working_model) {
  function(t, beta, X) {
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

#' @importFrom torch torch_reshape autograd_grad torch_cat
#' @importFrom purrr map
#' @noRd
ddL <- function(Lm, t, beta, X, p, weight = 1) {
  beta$grad <- NULL
  l <- dL(Lm, t, beta, X, weight = weight)[[1]]
  r <- purrr::map(1:length(l), function(index) {
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
    g <- torch::autograd_grad(
      l[index],
      list(t),
      retain_graph = TRUE,
      create_graph = TRUE,
      allow_unused = TRUE
    )[[1]]

    if(is.null(g)) {
      torch::torch_zeros(c(1, k))
    }
    else {
      torch::torch_reshape(g, shape = c(1, k))
    }
  })
  torch::torch_cat(r) # p x K
}

#' @importFrom torch optim_lbfgs torch_tensor
#' @noRd
B <- function(Lm_fn, psi, design_matrix, Q, p, init = NULL) {
  if(is.null(init)) init <- rep(0, p)
  beta <- torch_tensor(init, requires_grad = TRUE)
  optimizer <- torch::optim_lbfgs(beta, line_search_fn = "strong_wolfe")
  for (iter in 1:2) {
    optimizer$step(function() {
      optimizer$zero_grad()
      value <- Lm_fn(psi, beta, design_matrix)$mul(Q)$sum()
      value$backward()
      value
    })
    optimizer$zero_grad()
    if (max(abs(as.numeric(beta$grad))) < 1e-8) break
  }
  beta
}

#' @noRd
objective <- function(Lm, psi, beta, design_matrix, Q) {
  dL(Lm, psi, beta, design_matrix, weight = Q)[[1]]
}

#' @importFrom torch torch_cat autograd_grad
#' @importFrom purrr map
#' @noRd
dobjective_dpsi <- function(Lm, psi, beta, design_matrix, Q, p) {
  n <- dim(design_matrix)[1]
  k <- dim(design_matrix)[2]
  obj <- objective(Lm, psi, beta, design_matrix, Q)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], psi, retain_graph = TRUE)[[1]]$reshape(c(1, n, k))
  })
  return(torch::torch_cat(r)$transpose(1, 3)) # K x n x p
}

#' @importFrom torch autograd_grad torch_cat
#' @importFrom purrr map
#' @noRd
dobjective_dQ <- function(Lm, psi, beta, design_matrix, Q, p) {
  n <- dim(design_matrix)[1]
  obj <- objective(Lm, psi, beta, design_matrix, Q)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], Q, retain_graph = TRUE)[[1]]$reshape(c(n, 1))
  })
  torch::torch_cat(r, dim = 2) # n x p
}

#' @importFrom purrr map
#' @importFrom torch autograd_grad torch_cat
#' @noRd
dobjective_dbeta <- function(Lm, psi, beta, design_matrix, Q, p) {
  obj <- objective(Lm, psi, beta, design_matrix, Q)
  r <- map(1:p, function(index) {
    autograd_grad(obj[index], beta, retain_graph = TRUE)[[1]]$reshape(c(p, 1))
  })
  torch::torch_cat(r, dim = 2) # p x p
}

#' Invert E_Q\[Lddot\], or surface that it's degenerate
#'
#' Returns NULL when the objective Hessian is singular or ill-conditioned.
#' @noRd
invert_objective_hessian <- function(H, tol = 1e-10) {
  Hm <- as.matrix(H$detach())
  if(!all(is.finite(Hm))) return(NULL)

  rc <- tryCatch(rcond(Hm), error = function(e) NA_real_)
  if(!is.finite(rc) || rc < tol) return(NULL)
  inv <- tryCatch(solve(Hm), error = function(e) NULL)
  if(is.null(inv) || !all(is.finite(inv))) return(NULL)
  list(inv = torch::torch_tensor(inv), rcond = rc)
}

dB_dpsi <- function(Lm, psi, beta, design_matrix, Q, p, Hinv = NULL) {
  if(is.null(Hinv)) {
    check <- invert_objective_hessian(dobjective_dbeta(Lm, psi, beta, design_matrix, Q, p))
    if(is.null(check)) {
      stop("E_Q[Lddot] is singular; cannot form dB_dpsi.")
    }
    Hinv <- check$inv
  }
  -dobjective_dpsi(Lm, psi, beta, design_matrix, Q, p)$matmul(Hinv)
}

dB_dQ <- function(Lm, psi, beta, design_matrix, Q, p, Hinv = NULL) {
  if(is.null(Hinv)) {
    check <- invert_objective_hessian(dobjective_dbeta(Lm, psi, beta, design_matrix, Q, p))
    if(is.null(check)) {
      stop("E_Q[Lddot] is singular; cannot form dB_dpsi.")
    }
    Hinv <- check$inv
  }
  -dobjective_dQ(Lm, psi, beta, design_matrix, Q, p)$matmul(Hinv)
}

#' @importFrom torch torch_tensor
#' @noRd
normalizing_matrix <- function(Lm, psi, beta, design_matrix, Q, p) {
  n <- dim(design_matrix)[1]
  M <- ddL(Lm, psi, beta, design_matrix, p, Q)

  Mm <- as.matrix(M$detach())
  k <- kappa(Mm, exact = TRUE)

  if(!is.finite(k) || k > 1e12) {
    stop("Normalizing matrix M is numerically singular ",
         "(condition number ", signif(k, 3), "). ",
         "Check the working-model design for collinearity.", call. = FALSE)
  }
  if(k > 1e8) {
    warning("Normalizing matrix is ill-conditioned (condition number ", signif(k, 3), ").", call. = FALSE)
  }

  torch::torch_tensor(M$inverse())
}

#' Per-observation gradient of the loss wrt beta (the \dot{L}_m / D_2 term)
#' Returns an n x p tensor whose i-th row is dL_m/dbeta for observation i
#' @noRd
batched_dL_dbeta <- function(Lm, psi, beta, design_matrix, p, batched = TRUE) {
  n <- dim(design_matrix)[1]

  if(batched) {
    B <- beta$detach()$unsqueeze(1)$expand(c(n, p))$clone()$requires_grad_(TRUE)
    return(torch::autograd_grad(Lm(psi, B, design_matrix)$sum(), B, retain_graph = TRUE)[[1]])
  }

  Lm_vec <- Lm(psi, beta, design_matrix)
  torch::torch_cat(
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
eif <- function(Lm, psi, beta, design_matrix, Q, Delta, p, batched_beta) {
  n <- dim(design_matrix)[1]
  k <- dim(design_matrix)[2]
  Minv <- normalizing_matrix(Lm, psi, beta, design_matrix, Q, p)
  D2 <- batched_dL_dbeta(Lm, psi, beta, design_matrix, p, batched_beta)

  grad_blocks <- batched_NablaLdot(Lm, psi, beta, design_matrix, p, k)
  D1 <- torch::torch_cat(
    purrr::map(grad_blocks, function(gb) {
      (gb * Delta)$sum(dim = 2)$reshape(c(n, 1))
    }),
    dim = 2
  )

  eif <- as.matrix((D1 + D2)$matmul(-Minv))
  return(eif)
}
