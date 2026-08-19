test_that("out-of-support proposals aren't fatal", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)
  it <- cate_internals(data, nu, formula = ~X4, tmle = TRUE)

  tal <- new_bayes_tally()
  out <- bayes_log_density(
    torch::torch_tensor(rep(1e6, it$problem$p)),
    it$problem,
    it$spec,
    it$fit$state,
    it$fit$final$clever,
    it$fit$final$K_Q,
    it$condvar,
    prior = it$prior,
    control = list(eps_max = 25, min_ess = 2 * it$problem$p),
    tally = tal
  )

  expect_equal(out$log.density, -Inf)
  expect_length(out$beta, it$problem$p)
  expect_equal(unname(tally_counts(tal)["eps_max"]), 1L)
})

test_that("guards are active inside run_bayes_tmle", {
  data <- sim1_data(n = 200L)
  nu <- oracle_nuisance_cate(data)

  expect_warning(fit <- cate(data, X = paste0("X", 1:4), A = "A", Y = "Y", formula = ~X4,
              nuisance_estimates = nu,
              bayes = bayes_control(chains = 1L, draws = 100L, warmup = 100L,
                                    scale = 1e3, seed = 1L)))

  expect_true(sum(fit$bayes_tmle$rejected) > 0)
  expect_true(all(is.finite(as.array(fit$bayes_tmle$draws))))
})
