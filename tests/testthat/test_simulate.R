context("simulate")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# add covariate effects
ntime <- 35
xsim <- rnorm(ntime)
cov_fn <- function(popmat, x) {
  popmat[survival(popmat)] <- popmat[survival(popmat)] * plogis(x)
  popmat
}
cov_eff <- covariates(x = xsim,
                      fun = cov_fn)

# add density dependence
dd_masks <- list(reproduction(mat, dims = 4:5))
dd_fns <- list(
  function(x, n) x * (1 / (1 + 1e-4 * sum(n)))
)
dd <- density_dependence(masks = dd_masks, funs = dd_fns)

# add demographic stochasticity
demo_mask <- all_stages(mat)
demo_fn <- function(x) x + 1
demostoch <- demographic_stochasticity(masks = demo_mask, funs = demo_fn)

# add environmental stochasticity
mask_list <- list(reproduction(mat, dims = 4:5), survival(mat))
fn_list <- list(
  function(x) x + 2,
  function(x) x + 0.01
)
envstoch <- environmental_stochasticity(masks = mask_list, funs = fn_list)

# add rescaling method (alternative form of density dependence)
rescale_mask <- all_stages(mat)
rescale_fn <- function(x) 200 * (x / sum(x))
resc <- density_dependence_n(rescale_mask, rescale_fn)

# define simulation settings
nsim <- 10
ntime <- 50

test_that("simulate returns correct abundances without extra processes", {

  # simulate trajectories with no uncertainty
  dyn <- dynamics(mat)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(ntime))
    target[, , (i + 1)] <- floor(target[, , i] %*% t(dyn$matrix))
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # are lambdas roughly equal to leading eigenvalue?
  lambda_manual <-
    mean(apply(target[, , ntime], 1, sum) /
         apply(target[, , ntime - 1], 1, sum))
  lambda_value <-
    mean(apply(value[, , ntime], 1, sum) /
           apply(value[, , ntime - 1], 1, sum))
  lambda_real <- eigen(mat)$values[1]
  expect_equal(round(Re(lambda_real), 4),
               round(lambda_manual, 4))
  expect_equal(round(Re(lambda_real), 4),
               round(lambda_value, 4))

})

test_that("simulate returns correct abundances with covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(dyn$ntime))
    target[, , i + 1] <- floor(target[, , i] %*% t(dyn$matrix[[i]]))
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with density dependence and covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff, dd)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  dd_manual <- function(x, n) {
    x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
    x
  }
  for (i in seq_len(dyn$ntime)) {
    mat_tmp <- lapply(seq_len(nsim), function(x) t(dd_manual(dyn$matrix[[i]], target[x, , i])))
    target[, , i + 1] <- floor(t(mapply(`%*%`, lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp)))
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with demographic stochasticity and covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff, demostoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(dyn$ntime))
    target[, , i + 1] <- floor(target[, , i] %*% t(dyn$matrix[[i]])) + 1
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with demographic and environmental stochasticity and covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff, demostoch, envstoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(dyn$ntime)) {
    mat_tmp <- dyn$matrix[[i]]
    mat_tmp[1, 4:5] <- mat_tmp[1, 4:5] + 2
    idx <- row(mat_tmp) == col(mat_tmp)
    mat_tmp[idx] <- mat_tmp[idx] + 0.01
    target[, , i + 1] <- floor(target[, , i] %*% t(mat_tmp)) + 1
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with rescale density dependence and covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff, resc)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(dyn$ntime)) {
    target[, , i + 1] <- target[, , i] %*% t(dyn$matrix[[i]])
    target[, , i + 1] <- floor(t(apply(target[, , i + 1], 1, function(x) 200 * (x / sum(x)))))
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with density dependence", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, dd)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  dd_manual <- function(x, n) {
    x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
    x
  }
  for (i in seq_len(ntime)) {
    mat_tmp <- lapply(seq_len(nsim), function(x) t(dd_manual(dyn$matrix, target[x, , i])))
    target[, , i + 1] <- floor(t(mapply(`%*%`, lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp)))
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with density dependence and demostoch", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, dd, demostoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  dd_manual <- function(x, n) {
    x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
    x
  }
  for (i in seq_len(ntime)) {
    mat_tmp <- lapply(seq_len(nsim), function(x) t(dd_manual(dyn$matrix, target[x, , i])))
    target[, , i + 1] <- floor(t(mapply(`%*%`, lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp)) + 1)
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with density dependence and envstoch", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, dd, envstoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  dd_manual <- function(x, n) {
    x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
    x
  }
  for (i in seq_len(ntime)) {
    mat_tmp <- dyn$matrix
    mat_tmp[1, 4:5] <- mat_tmp[1, 4:5] + 2
    idx <- row(mat_tmp) == col(mat_tmp)
    mat_tmp[idx] <- mat_tmp[idx] + 0.01
    mat_tmp <- lapply(seq_len(nsim), function(x) t(dd_manual(mat_tmp, target[x, , i])))
    target[, , i + 1] <- floor(t(mapply(`%*%`, lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp)))
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})


test_that("simulate errors informatively when initial values have unsuitable dims", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, dd, envstoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  expect_error(
    simulate(
      dyn,
      nsim = nsim,
      init = init_set[-1],
      options = list(ntime = ntime, tidy_abundances = floor)
    ),
    "must have dimensions \\(10,5\\) or \\(5\\)"
  )

})
