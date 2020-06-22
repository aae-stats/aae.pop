context("dynamics")

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
demo_mask <- all_classes(mat)
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
rescale_mask <- all_classes(mat)
rescale_fn <- function(x) 200 * (x / sum(x))
resc <- density_dependence_n(rescale_mask, rescale_fn)

test_that("dynamics object returns correct processes for different combinations of processes", {

  # test basic dynamics object with no additional processes
  dyn_obj <- dynamics(mat)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, 1L)
  expect_null(dyn_obj$covariates)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_null(dyn_obj$density_dependence)
  expect_null(dyn_obj$density_dependence_n)

  # test basic dynamics object with covariates
  dyn_obj <- dynamics(mat, cov_eff)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, cov_eff$ntime)
  expect_equal(dyn_obj$covariates, cov_eff)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_null(dyn_obj$density_dependence)
  expect_null(dyn_obj$density_dependence_n)

  # test basic dynamics object with environmental stochasticity
  dyn_obj <- dynamics(mat, envstoch)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, 1L)
  expect_null(dyn_obj$covariates)
  expect_equal(dyn_obj$environmental_stochasticity, envstoch)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_null(dyn_obj$density_dependence)
  expect_null(dyn_obj$density_dependence_n)

  # test basic dynamics object with demographic stochasticity
  dyn_obj <- dynamics(mat, demostoch)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, 1L)
  expect_null(dyn_obj$covariates)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_equal(dyn_obj$demographic_stochasticity, demostoch)
  expect_null(dyn_obj$density_dependence)
  expect_null(dyn_obj$density_dependence_n)

  # test basic dynamics object with density_dependence
  dyn_obj <- dynamics(mat, dd)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, 1L)
  expect_null(dyn_obj$covariates)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_equal(dyn_obj$density_dependence, dd)
  expect_null(dyn_obj$density_dependence_n)

  # test basic dynamics object with rescale density dependence
  dyn_obj <- dynamics(mat, resc)
  expect_equal(dyn_obj$matrix, mat)
  expect_equal(dyn_obj$nclass, ncol(mat))
  expect_equal(dyn_obj$nspecies, 1L)
  expect_equal(dyn_obj$ntime, 1L)
  expect_null(dyn_obj$covariates)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_null(dyn_obj$density_dependence)
  expect_equal(dyn_obj$density_dependence_n, resc)

})


test_that("dynamics objects can be updated", {

  # create basic object
  dyn_obj <- dynamics(mat)

  # update by adding covariates
  dyn_new <- update(dyn_obj, cov_eff)
  expect_null(dyn_obj$covariates)
  expect_equal(dyn_new$covariates, cov_eff)

  # update by adding density dependence
  dyn_new <- update(dyn_obj, dd)
  expect_null(dyn_obj$density_dependence)
  expect_equal(dyn_new$density_dependence, dd)

  # update by adding demographic stochasticity
  dyn_new <- update(dyn_obj, demostoch)
  expect_null(dyn_obj$demographic_stochasticity)
  expect_equal(dyn_new$demographic_stochasticity, demostoch)

  # update by adding environmental stochasticity
  dyn_new <- update(dyn_obj, envstoch)
  expect_null(dyn_obj$environmental_stochasticity)
  expect_equal(dyn_new$environmental_stochasticity, envstoch)

  # update by adding rescale density depdencen
  dyn_new <- update(dyn_obj, resc)
  expect_null(dyn_obj$density_dependence_n)
  expect_equal(dyn_new$density_dependence_n, resc)

  # update by changing covariates
  dyn_new <- update(dyn_obj, cov_eff)
  expect_equivalent(dyn_new$covariates$x, xsim)
  xnew <- rnorm(length(xsim))
  dyn_new <- update(dyn_new, covariates(x = xnew, fun = cov_fn))
  expect_equivalent(dyn_new$covariates$x, xnew)

})

test_that("dynamics object errors informatively with unsuitable processes", {

  # no matrix passed to dynamics
  expect_error(dynamics(),
               "matrix must be provided")

  # a non-matrix passed as first argument
  expect_error(dynamics(dd),
               "matrix must be a two-dimensional array or matrix")

  # non-process passed to dynamics
  expect_error(
    dynamics(mat, rnorm(10)),
    "must be one of covariates, environmental_stochasticity"
  )

  # multiples of the same process
  expect_error(
    dynamics(mat, dd, dd),
    "Multiple objects provided for the following processes: density_dependence"
  )

  # multiples of the same two processes
  expect_error(
    dynamics(mat, dd, dd, demostoch, demostoch),
    "Multiple objects provided for the following processes: demographic_stochasticity, and density_dependence"
  )

})

test_that("dynamics objects can be plotted", {

  # create basic object
  dyn_obj <- dynamics(mat)

  # plot it
  expect_silent(plot(dyn_obj))

})
