context("simulate")

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
xsim2 <- rnorm(ntime)
cov_fn <- function(mat, x) {
  mat * plogis(x)
}
cov_eff <- covariates(
  masks = survival(mat),
  funs = cov_fn
)

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

# add extra versions of the above with either lists or no lists (to
#   catch the things not already tested above)
cov_eff2 <- covariates(
  masks = list(survival(mat, dims = 3:5), survival(mat, dims = 1)),
  funs = list(cov_fn, cov_fn)
)
rep_cov_eff2 <- replicated_covariates(
  masks = list(survival(mat, dims = 3:5), survival(mat, dims = 1)),
  funs = list(cov_fn, cov_fn)
)
dd2 <- density_dependence(
  masks = reproduction(mat, dims = 4:5),
  funs = \(x, n) x * (1 / (1 + 1e-4 * sum(n)))
)
demostoch2 <- demographic_stochasticity(
  masks = list(demo_mask, demo_mask),
  funs = list(demo_fn, demo_fn)
)
envstoch2 <- environmental_stochasticity(
  masks = reproduction(mat, dims = 4:5),
  funs = \(x) x + 2
)
resc2 <- density_dependence_n(
  list(rescale_mask, rescale_mask),
  list(rescale_fn, rescale_fn)
)
resc_pre2 <- add_remove_pre(
  list(rescale_mask, rescale_mask),
  list(rescale_fn, rescale_fn)
)
resc_post2 <- add_remove_post(
  list(rescale_mask, rescale_mask),
  list(rescale_fn, rescale_fn)
)

# extra covariates term to test aux and named vars in format_covariates
cov_eff3 <- covariates(
  masks = survival(mat),
  funs = \(mat, x, y) mat * plogis(x + y)
)
cov_eff4 <- covariates(
  masks = survival(mat),
  funs = \(mat, x) mat * plogis(x["x"] + x["y"])
)

# set pre-cooked dd functions for testing
dd_bh <- density_dependence(
  masks = reproduction(mat, dims = 4:5),
  funs = beverton_holt(k = 100)
)
dd_ricker <- density_dependence(
  masks = reproduction(mat, dims = 4:5),
  funs = ricker(k = 100)
)

# define default imulation settings
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
  expect_equal(round(Re(lambda_real), 2),
               round(lambda_manual, 2))
  expect_equal(round(Re(lambda_real), 2),
               round(lambda_value, 2))

})

test_that("simulate returns correct abundances with covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(
      covariates = list(
        lapply(seq_along(xsim), function(i) xsim[i])
      )
    )
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_along(xsim)) {
    target[, , i + 1] <- floor(
      target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
    )
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # test again with format_covariates
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(covariates = format_covariates(xsim))
  )
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with named covariates
          and aux vars", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff3)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(
                covariates = format_covariates(xsim, names = "x", aux = 1)
              )
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <- floor(
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i], 1))
              )
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with named covariates
          in mixed orders", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff4)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(
                covariates = format_covariates(
                  cbind(xsim, xsim2), names = c("x", "y")
                )
              )
            )
            value2 <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(
                covariates = format_covariates(
                  cbind(xsim2, xsim), names = c("y", "x")
                )
              )
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <- floor(
                target[, , i] %*% t(
                  dyn$covariates(dyn$matrix, c(x = xsim[i], y = xsim2[i]))
                )
              )
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)
            expect_equal(target, value2)

          })

test_that("simulate returns correct abundances with list covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, cov_eff2)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(
      covariates = list(
        lapply(seq_along(xsim), function(i) xsim[i])
      )
    )
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_along(xsim)) {
    target[, , i + 1] <- floor(
      target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
    )
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # test again with format_covariates
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(covariates = format_covariates(xsim))
  )
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with replicated_covariates", {

  # simulate trajectories with no uncertainty but with covariates
  dyn <- dynamics(mat, rep_cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(
      replicated_covariates = list(
        lapply(seq_len(nrow(xsim_rep)), function(i) xsim_rep[i, ])
      )
    )
  )
  target <- array(dim = dim(value))
  target[, , 1] <- init_set
  for (i in seq_len(nrow(xsim_rep))) {
    for (j in seq_len(ncol(xsim_rep))) {
      target[j, , i + 1] <- floor(
        target[j, , i] %*%
          t(dyn$replicated_covariates(dyn$matrix, xsim_rep[i, j]))
      )
    }
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

  # test again with format_covariates
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(replicated_covariates = format_covariates(xsim_rep))
  )
  expect_equal(target, value)

})

test_that(
          "simulate returns correct abundances with list replicated_covariates",
          {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, rep_cov_eff2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(
                replicated_covariates = list(
                  lapply(seq_len(nrow(xsim_rep)), function(i) xsim_rep[i, ])
                )
              )
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_len(nrow(xsim_rep))) {
              for (j in seq_len(ncol(xsim_rep))) {
                target[j, , i + 1] <- floor(
                  target[j, , i] %*%
                    t(dyn$replicated_covariates(dyn$matrix, xsim_rep[i, j]))
                )
              }
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

            # test again with format_covariates
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(replicated_covariates = format_covariates(xsim_rep))
            )
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances
           with density dependence and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, dd)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            dd_manual <- function(x, n) {
              x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
              x
            }
            for (i in seq_along(xsim)) {
              mat_tmp <- lapply(
                seq_len(nsim),
                function(x) {
                  t(dd_manual(
                    dyn$covariates(dyn$matrix, xsim[i]),
                    target[x, , i]
                  ))
                }
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp
                )
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances
           with list density dependence and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, dd2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            dd_manual <- function(x, n) {
              x[1, 4:5] <- dd_fns[[1]](x[1, 4:5], n)
              x
            }
            for (i in seq_along(xsim)) {
              mat_tmp <- lapply(
                seq_len(nsim),
                function(x) {
                  t(dd_manual(
                    dyn$covariates(dyn$matrix, xsim[i]), target[x, , i]
                  ))
                }
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp
                )
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           demographic stochasticity and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, demostoch)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <- floor(
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              ) + 1
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           list demographic stochasticity and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, demostoch2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <- floor(
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              ) + 2
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with demographic
           and environmental stochasticity and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, demostoch, envstoch)
            init_set <- matrix(
              rpois(nstage * nsim, lambda = 20), ncol = nstage
            )
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              mat_tmp <- dyn$covariates(dyn$matrix, xsim[i])
              mat_tmp[1, 4:5] <- mat_tmp[1, 4:5] + 2
              idx <- row(mat_tmp) == col(mat_tmp)
              mat_tmp[idx] <- mat_tmp[idx] + 0.01
              target[, , i + 1] <- floor(target[, , i] %*% t(mat_tmp)) + 1
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with demographic
           and list environmental stochasticity and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, demostoch, envstoch2)
            init_set <- matrix(
              rpois(nstage * nsim, lambda = 20), ncol = nstage
            )
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              mat_tmp <- dyn$covariates(dyn$matrix, xsim[i])
              mat_tmp[1, 4:5] <- mat_tmp[1, 4:5] + 2
              target[, , i + 1] <- floor(target[, , i] %*% t(mat_tmp)) + 1
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           rescale density dependence and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc)
            init_set <- matrix(
              rpois(nstage * nsim, lambda = 20), ncol = nstage
            )
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <-
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(t(
                apply(target[, , i + 1], 1, function(x) 200 * (x / sum(x)))
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           list rescale density dependence and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <-
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(t(
                apply(target[, , i + 1], 1, function(x) 200 * (x / sum(x)))
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           removals pre-update and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc_pre)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              tmp_target <- t(
                apply(target[, , i], 1, function(x) 200 * (x / sum(x)))
              )
              target[, , i + 1] <-
                tmp_target %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(target[, , i + 1])
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           list removals pre-update and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc_pre2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              tmp_target <- t(
                apply(target[, , i], 1, function(x) 200 * (x / sum(x)))
              )
              target[, , i + 1] <-
                tmp_target %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(target[, , i + 1])
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           removals post-update and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc_post)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <-
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(t(
                apply(target[, , i + 1], 1, function(x) 200 * (x / sum(x)))
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           list removals post-update and covariates", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, cov_eff, resc_post2)
            init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
            value <- simulate(
              dyn,
              nsim = nsim,
              init = init_set,
              options = list(ntime = ntime, tidy_abundances = floor),
              args = list(covariates = format_covariates(xsim))
            )
            target <- array(dim = dim(value))
            target[, , 1] <- init_set
            for (i in seq_along(xsim)) {
              target[, , i + 1] <-
                target[, , i] %*% t(dyn$covariates(dyn$matrix, xsim[i]))
              target[, , i + 1] <- floor(t(
                apply(target[, , i + 1], 1, function(x) 200 * (x / sum(x)))
              ))
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
    mat_tmp <- lapply(
      seq_len(nsim), function(x) t(dd_manual(dyn$matrix, target[x, , i]))
    )
    target[, , i + 1] <- floor(t(
      mapply(
        `%*%`, lapply(seq_len(nsim), function(x) target[x, , i]), mat_tmp
      )
    ))
  }
  class(target) <- c("simulation", "array")
  expect_equal(target, value)

})

test_that("simulate returns correct abundances with
           density dependence and demostoch", {

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
              mat_tmp <- lapply(
                seq_len(nsim),
                \(x) t(dd_manual(dyn$matrix, target[x, , i]))
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]),
                  mat_tmp
                )
              ) + 1)
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           density dependence and envstoch", {

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
              mat_tmp <- lapply(
                seq_len(nsim), function(x) t(dd_manual(mat_tmp, target[x, , i]))
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]),
                  mat_tmp
                )
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           bh density dependence", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, dd_bh)
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
              x[1, 4:5] <- x[1, 4:5] / (1 + x[1, 4:5] * sum(n) / 100)
              x
            }
            for (i in seq_len(ntime)) {
              mat_tmp <- dyn$matrix
              mat_tmp <- lapply(
                seq_len(nsim), function(x) t(dd_manual(mat_tmp, target[x, , i]))
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]),
                  mat_tmp
                )
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns correct abundances with
           ricker density dependence", {

            # simulate trajectories with no uncertainty but with covariates
            dyn <- dynamics(mat, dd_ricker)
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
              x[1, 4:5] <- x[1, 4:5] * exp(1 - sum(n) / 100) / exp(1)
              x
            }
            for (i in seq_len(ntime)) {
              mat_tmp <- dyn$matrix
              mat_tmp <- lapply(
                seq_len(nsim), function(x) t(dd_manual(mat_tmp, target[x, , i]))
              )
              target[, , i + 1] <- floor(t(
                mapply(
                  `%*%`,
                  lapply(seq_len(nsim), function(x) target[x, , i]),
                  mat_tmp
                )
              ))
            }
            class(target) <- c("simulation", "array")
            expect_equal(target, value)

          })

test_that("simulate returns reproducible outputs when seed is set", {

  # update environmental stochasticity to actually be stochastic
  mask_list_stoch <- list(reproduction(mat, dims = 4:5), survival(mat))
  fn_list_stoch <- list(
    function(x) rpois(length(x), lambda = x),
    function(x) rnorm(length(x), mean = x, sd = 0.01)
  )
  envstoch_stoch <- environmental_stochasticity(
    masks = mask_list_stoch, funs = fn_list_stoch
  )

  # simulate two trajectories with noise but same seed
  dyn <- dynamics(mat, cov_eff, dd, envstoch_stoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    seed = 123,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(covariates = format_covariates(xsim))
  )
  target <- simulate(
    dyn,
    nsim = nsim,
    seed = 123,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(covariates = format_covariates(xsim))
  )
  expect_equal(target, value)

  # and now with a different seed
  target <- simulate(
    dyn,
    nsim = nsim,
    seed = 124,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor),
    args = list(covariates = format_covariates(xsim))
  )
  expect_false(all(target == value))

})

test_that("simulate errors informatively when inits have unsuitable dims", {

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

test_that("simulate errors informatively when dyn args have unsuitable dims", {

  # when multiple dynamic arguments are passed with conflicting dimensions
  dyn <- dynamics(mat, cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  expect_error(
    simulate(
      dyn,
      nsim = nsim,
      init = init_set,
      options = list(ntime = ntime, tidy_abundances = floor),
      args = list(
        covariates = list(
          lapply(seq_len(length(xsim) - 1L), function(i) xsim[i]),
          lapply(seq_along(xsim), function(i) xsim[i])
        )
      )
    ),
    "must have the same length"
  )

  # when replicated arguments have the wrong dimensions
  dyn <- dynamics(mat, rep_cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  expect_error(
    simulate(
      dyn,
      nsim = nsim,
      init = init_set,
      options = list(ntime = ntime, tidy_abundances = floor),
      args = list(
        replicated_covariates = list(
          lapply(seq_len(nrow(xsim_rep)), function(i) xsim_rep[i, 1:5])
        )
      )
    ),
    "should have nsim columns"
  )

  # when dynamic and replicated arguments have mismatched dimensions
  dyn <- dynamics(mat, cov_eff, rep_cov_eff)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  expect_error(
    simulate(
      dyn,
      nsim = nsim,
      init = init_set,
      options = list(ntime = ntime, tidy_abundances = floor),
      args = list(
        covariates = list(
          lapply(seq_len(length(xsim) - 1L), function(i) xsim[i])
        ),
        replicated_covariates = list(
          lapply(seq_len(nrow(xsim_rep)), function(i) xsim_rep[i, ])
        )
      )
    ),
    "must have the same length"
  )

})

test_that("simulation objects can be plotted and printed", {

  # basic simulation
  dyn <- dynamics(mat, dd, demostoch)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )

  # plot it
  expect_silent(plot(value))

  # and check other outputs
  expect_output(summary(value), "Simulated population has a")
  expect_output(print(value), "Simulated population dynamics for a")
  expect_true(is.simulation(value))
  expect_false(is.simulation(1L))

})
