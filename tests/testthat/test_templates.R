context("templates")

test_that("templates return working dynamics objects", {

  # simulate from a Murray cod object
  dyn <- murraycod()
  sim <- simulate(dyn)
  expect_equal(dim(sim), c(1L, 25L, 51L))

  # expect most processed to be defined
  expect_null(dyn$density_dependence_n)
  expect_equal(class(dyn$covariates), c("covariates", "function"))
  expect_equal(class(dyn$density_dependence), c("density_dependence", "function"))
  expect_equal(class(dyn$demographic_stochasticity), c("demographic_stochasticity", "function"))
  expect_equal(class(dyn$environmental_stochasticity), c("environmental_stochasticity", "function"))

})
