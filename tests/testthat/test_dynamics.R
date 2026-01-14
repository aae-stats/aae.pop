context("dynamics")

# setup: simulate some data to test with
nsim <- 10
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# add covariate effects
ntime <- 35
xsim <- rnorm(ntime)
cov_fn <- function(mat, x) {
  mat * plogis(x)
}
cov_eff <- covariates(masks = survival(mat),
                      funs = cov_fn)

# add replicated_covariates effects
xsim_rep <- matrix(rnorm(ntime * nsim), ncol = nsim)
rep_cov_eff <- replicated_covariates(
  masks = survival(mat),
  funs = cov_fn
)

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
resc_pre <- add_remove_pre(rescale_mask, rescale_fn)
resc_post <- add_remove_post(rescale_mask, rescale_fn)

test_that("dynamics object returns correct processes
           for different combinations of processes", {

             # test basic dynamics object with no additional processes
             dyn_obj <- dynamics(mat)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with covariates
             dyn_obj <- dynamics(mat, cov_eff)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_equal(dyn_obj$covariates, cov_eff)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with replicated_covariates
             dyn_obj <- dynamics(mat, rep_cov_eff)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_equal(dyn_obj$replicated_covariates, rep_cov_eff)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with environmental stochasticity
             dyn_obj <- dynamics(mat, envstoch)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_equal(dyn_obj$environmental_stochasticity, envstoch)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with demographic stochasticity
             dyn_obj <- dynamics(mat, demostoch)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_equal(dyn_obj$demographic_stochasticity, demostoch)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with density_dependence
             dyn_obj <- dynamics(mat, dd)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_equal(dyn_obj$density_dependence, dd)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with rescale density dependence
             dyn_obj <- dynamics(mat, resc)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_equal(dyn_obj$density_dependence_n, resc)

             # test basic dynamics object with removals pre-update
             dyn_obj <- dynamics(mat, resc_pre)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_equal(dyn_obj$add_remove_pre, resc_pre)
             expect_null(dyn_obj$add_remove_post)
             expect_null(dyn_obj$density_dependence_n)

             # test basic dynamics object with removals post-update
             dyn_obj <- dynamics(mat, resc_post)
             expect_equal(dyn_obj$matrix, mat)
             expect_equal(dyn_obj$nclass, ncol(mat))
             expect_null(dyn_obj$covariates)
             expect_null(dyn_obj$replicated_covariates)
             expect_null(dyn_obj$environmental_stochasticity)
             expect_null(dyn_obj$demographic_stochasticity)
             expect_null(dyn_obj$density_dependence)
             expect_null(dyn_obj$add_remove_pre)
             expect_equal(dyn_obj$add_remove_post, resc_post)
             expect_null(dyn_obj$density_dependence_n)

           })


test_that("dynamics objects can be updated", {

  # create basic object
  dyn_obj <- dynamics(mat)

  # update by adding covariates
  dyn_new <- update(dyn_obj, cov_eff)
  expect_equal(class(dyn_new$covariates), c("covariates", "function"))
  expect_equal(dyn_new$covariates, cov_eff)

  # update by adding replicated_covariates
  dyn_new <- update(dyn_obj, rep_cov_eff)
  expect_equal(
    class(dyn_new$replicated_covariates),
    c("replicated_covariates", "function")
  )
  expect_equal(dyn_new$replicated_covariates, rep_cov_eff)

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

  # update by changing covariate function
  dyn_new <- update(dyn_obj, cov_eff)
  new_mat <- dyn_new$matrix
  new_mat[survival(new_mat)] <- plogis(10) * new_mat[survival(new_mat)]
  expect_equal(dyn_new$covariates(dyn_new$matrix, 10), new_mat)
  new_fn <- function(mat, x) mat + x
  dyn_new <- update(
    dyn_new,
    covariates(masks = survival(dyn_new$matrix), funs = new_fn)
  )
  new_mat <- dyn_new$matrix
  new_mat[survival(new_mat)] <- 0.01 + new_mat[survival(new_mat)]
  expect_equal(dyn_new$covariates(dyn_new$matrix, 0.01), new_mat)

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
    "must be one of covariates, replicated_covariates"
  )

  # multiples of the same process
  expect_error(
    dynamics(mat, dd, dd),
    "Multiple objects provided for the following processes: density_dependence"
  )

  # multiples of the same two processes
  expect_error(
    dynamics(mat, dd, dd, demostoch, demostoch),
    paste0("Multiple objects provided for the following processes:",
           " demographic_stochasticity, and density_dependence")
  )

})

test_that("processes error informatively with mismatched mask fn dims", {

  # no matrix passed to dynamics
  expect_error(
    density_dependence(
      masks = list(reproduction(mat)),
      funs = list(\(x) x + 1, \(x) x + 2)
    ),
    "must be one element of funs for each element of mask"
  )

})

test_that("dynamics objects can be plotted", {

  # create basic object
  dyn_obj <- dynamics(mat)

  # plot it
  expect_silent(plot(dyn_obj))

})
