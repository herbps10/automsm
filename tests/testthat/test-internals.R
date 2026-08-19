source(test_path("helper-invariants.R"))
source(test_path("helper-internals.R"))

test_that("internals helpers agree with exported wrappers", {
  data <- sim1_data(n = 150L)
  nu <- oracle_nuisance_cate(data)

  exported <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y", formula = ~X4,
                   tmle = tmle_control(), bayes = FALSE, nuisance_estimates = nu)
  it <- cate_internals(data, nu, formula = ~X4, tmle = TRUE, joint_seed = 1L)

  expect_equal(it$plugin, exported$plugin$est, tolerance = 0)
  expect_equal(it$base$onestep$est, exported$onestep$est, tolerance = 0)
  expect_equal(it$tmle$est, exported$tmle$est, tolerance = 0)
  expect_equal(it$problem$p, exported$p)
  expect_equal(it$problem$d, exported$d)
  expect_equal(it$problem$terms, exported$terms)
})

