test_that("simulating cate example data suns", {
  dat <- simulate_cate(N = 10, sigma = 0.1, seed = 1)
  expect_equal(names(dat), c("X1", "X2", "A", "Y"))
  expect_equal(nrow(dat), 10)
})

test_that("simulating dose-response example data runs", {
  dat <- simulate_dose_response(N = 100, treatments = 5, sigma = 0.1, seed = 1)
  expect_equal(names(dat), c("X1", "X2", "A", "Y"))
  expect_true(all(dat$A %in% 1:5))
  expect_equal(nrow(dat), 100)
})

test_that("simulating longitudinal dose-response example data runs", {
  dat <- simulate_longitudinal_dose_response(N = 500, tau = 2, seed = 1)
  expect_equal(names(dat), c("L1", "L2", "A1", "A2", "Y"))
  expect_equal(nrow(dat), 500)
})
