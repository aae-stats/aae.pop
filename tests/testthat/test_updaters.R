context("updaters")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))
dyn <- dynamics(mat)

nsim <- 10
ntime <- 20

test_that("updaters give correct outputs", {

  # test standard cross product
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  sim <- simulate(dyn, options = list(update = update_crossprod))
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime,
                   tidy_abundances = identity,
                   update = update_crossprod)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(ntime))
    target[, , (i + 1)] <- target[, , i] %*% t(dyn$matrix)
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # test the same with floored abundances
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime,
                   tidy_abundances = floor,
                   update = update_crossprod)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(ntime))
    target[, , (i + 1)] <- floor(target[, , i] %*% t(dyn$matrix))
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # test multinomial update
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime,
                   tidy_abundances = floor,
                   update = update_multinomial)
  )
  expect_equal(class(value), c("simulate", "array"))

  # test binomial update
  mat[survival(mat)] <- 0
  dyn <- dynamics(mat)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime,
                   tidy_abundances = floor,
                   update = update_binomial_leslie)
  )
  expect_equal(class(value), c("simulate", "array"))

})


test_that("updaters error informatively", {


})


