test_that("glm_gradient matches central differences", {
  for(linear in c(TRUE, FALSE)) {
    family <- fluctuation_family(linear)
    set.seed(1)

    g <- list(X = matrix(rnorm(300), 100, 3),
              offset = if(linear) rnorm(100) else stats::qlogis(runif(100, 0.1, 0.9)),
              target = if(linear) rnorm(100) else runif(100),
              weights = runif(100, 0.5, 3))
    for(qp in list(NULL, q_penalty(matrix(rnorm(300), 100, 3), rep(1/100, 100), 100))) {
      for(e in list(rep(0, 3), c(0.3, -0.2, 0.1), c(-1.5, 2.0, 0.4))) {
        expect_equal(
          glm_gradient(g, family, e, qp),
          numeric_jacobian(function(x) glm_objective(g, family, x, qp), e)[1, ],
          tolerance = 1e-6, info = paste("linear:", linear, "qpen:", !is.null(qp))
        )
      }
    }
  }
})

