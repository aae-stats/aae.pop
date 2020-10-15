context("rng")
set.seed(1242)

test_that("rng functions work as expected", {

  # only one function to test at this stage
  n <- 10
  mean_vals <- runif(n, min = 0.3, max = 0.7)
  sd_vals <- 0.1 * mean_vals

  # test without correlations
  value <- rmultiunit(n = 10000, mean = mean_vals, sd = sd_vals)
  mean_test <- apply(value, 2, mean)
  sd_test <- apply(value, 2, sd)
  expect_true(all(abs(mean_test - mean_vals) < 1e-2))
  expect_true(all(abs(sd_test - sd_vals) < 1e-2))

  # test with perfect correlations
  value <- rmultiunit(n = 10000, mean = mean_vals, sd = sd_vals, perfect_correlation = TRUE)
  mean_test <- apply(value, 2, mean)
  sd_test <- apply(value, 2, sd)
  expect_true(all(abs(mean_test - mean_vals) < 1e-2))
  expect_true(all(abs(sd_test - sd_vals) < 1e-2))
  expect_true(all(cor(value) > 0.99))

  # test with imperfect correlations
  sigma <- runif(n)
  sigma <- sigma %o% sigma
  sigma <- sigma + diag(x = 5, nrow = n, ncol = n)
  omega <- cov2cor(sigma)
  value <- rmultiunit(n = 100000, mean = mean_vals, sd = sd_vals, Sigma = sigma)
  mean_test <- apply(value, 2, mean)
  sd_test <- apply(value, 2, sd)
  expect_true(all(abs(mean_test - mean_vals) < 1e-2))
  expect_true(all(abs(sd_test - sd_vals) < 1e-2))
  expect_true(all(abs(cor(value) - omega) < 1e-2))

  # same but based on omega
  value <- rmultiunit(n = 100000, mean = mean_vals, sd = sd_vals, Omega = omega)
  mean_test <- apply(value, 2, mean)
  sd_test <- apply(value, 2, sd)
  expect_true(all(abs(mean_test - mean_vals) < 1e-2))
  expect_true(all(abs(sd_test - sd_vals) < 1e-2))
  expect_true(all(abs(cor(value) - omega) < 1e-2))

  # test behaviour if conflicting arguments provided
  value <- rmultiunit(n = 100000, mean = mean_vals, sd = sd_vals, Omega = omega, perfect_correlation = TRUE)
  mean_test <- apply(value, 2, mean)
  sd_test <- apply(value, 2, sd)
  expect_true(all(abs(mean_test - mean_vals) < 1e-2))
  expect_true(all(abs(sd_test - sd_vals) < 1e-2))
  expect_true(all(abs(cor(value) - omega) < 1e-2))

  # test warnings and errors
  expect_warning(rmultiunit(1, c(0.5, 0.5, 0.5), c(0.2, 0.2)),
                 "is not a multiple")
  expect_warning(rmultiunit(10, mean = c(0.5, 0.5), sd = c(0.2, 0.2, 0.2)),
                 "is not a multiple")
  expect_error(rmultiunit(n = 10, mean = rep(0.5, 3), sd = rep(0.2, 3), Sigma = sigma),
               "subscript out of bounds")

})
