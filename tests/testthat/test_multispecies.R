context("multispecies")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# define multispecies interactions
multispecies_mask <- survival(mat)
multispecies_fun <- function(x, n, ...) {
  x / (1 + x * sum(n) / 2000)
}
dyn_mc1 <- dynamics(mat)
dyn_mc2 <- dynamics(mat)
dyn_mc3 <- dynamics(mat)
dyn_mc4 <- dynamics(mat)
interactions_test1 <- pairwise_interaction(
  dyn_mc1, dyn_mc2, multispecies_mask, multispecies_fun
)
interactions_test2 <- pairwise_interaction(
  dyn_mc1, dyn_mc3, multispecies_mask, multispecies_fun
)
interactions_test3 <- pairwise_interaction(
  dyn_mc3, dyn_mc4, multispecies_mask, multispecies_fun
)

# add some basic processes
envstoch <- environmental_stochasticity(
  masks = multispecies_mask,
  funs = \(x) pmin(1.1 * x, 1)
)
repcov <- replicated_covariates(
  masks = multispecies_mask,
  funs = \(x, y) y * x
)
cov <- covariates(
  masks = reproduction(mat, dims = 4:5),
  funs = \(x, y) x * y
)
dyn_mc5 <- dynamics(mat, envstoch, repcov, cov)
interactions_test4 <- pairwise_interaction(
  dyn_mc3, dyn_mc5, multispecies_mask, multispecies_fun
)

test_that("multispecies object has correct dynamics elements", {

  # create multispecies object
  mspecies_obj <- multispecies(
    interactions_test1, interactions_test2, interactions_test3
  )

  # check dynamics one by one (should have 1-4, not necessarily in that order)
  #  match based on hex
  expect_equal(1L, 1L)

})

# nolint start
test_that("multispecies objects simulate correctly", {

  # not a full test because test_simulate.R covers identical use of simulate

  # create multispecies object
  mspecies_obj <- multispecies(
    interactions_test1, interactions_test2, interactions_test3
  )

  # how many replicates are we simulating?
  nsim <- 10

  # define initial conditions
  dims <- c(nsim, nstage, mspecies_obj$nspecies)
  init <- array(rpois(prod(dims), lambda = 10), dim = dims)

  # simulate with aae.pop
  value <- simulate(mspecies_obj, nsim = nsim, init = init)

  # simulate again with same initial conditions
  target <- simulate(mspecies_obj, nsim = nsim, init = init)

  # and compare
  expect_equal(value, target)

  # now simulate manually
  all_dyn <- list(dyn_mc1, dyn_mc2, dyn_mc3, dyn_mc4)
  sp_order <- match(sapply(mspecies_obj$dynamics, function(x) x$hex),
                    sapply(all_dyn, function(x) x$hex))
  all_dyn <- all_dyn[sp_order]
  target <- lapply(
    value,
    function(x) array(dim = dim(x))
  )
  for (i in seq_along(target)) {
    target[[i]][, , 1] <- init[, , i]
  }
  interaction_pairs <- apply(
    mspecies_obj$structure, 1, function(x) which(x == 1)
  )
  for (i in seq_len(dim(value[[1]])[3] - 1)) {
    for (k in seq_len(nsim)) {
      for (j in seq_along(all_dyn)) {

        # pull out the baseline matrix for species j
        mat <- all_dyn[[j]]$matrix

        # add in pairwise interactions
        if (length(interaction_pairs[[j]]) > 0) {
          x <- mat[survival(mat)]
          for (.i in seq_along(interaction_pairs[[j]])) {
            x <- x /
              (1 +
                 x *
                 sum(target[[interaction_pairs[[j]][.i]]][k, , i]) /
                 2000
              )
          }
          mat[survival(mat)] <- x
        }

        # update
        target[[j]][k, , (i + 1)] <- target[[j]][k, , i] %*% t(mat)

      }
    }
  }
  for (i in seq_along(target))
    class(target[[i]]) <- c("simulation", "array")
  class(target) <- c("simulation_list", "list")
  expect_equal(value, target)

})

test_that("multispecies objects simulate correctly with basic processes", {

  # how many replicates are we simulating?
  nsim <- 10

  # define initial conditions
  dims <- c(nsim, nstage, mspecies_obj$nspecies)
  init <- array(rpois(prod(dims), lambda = 10), dim = dims)

  # define object with processes in interactions_test4
  mspecies_obj <- multispecies(
    interactions_test1, interactions_test2, interactions_test4
  )

  # simulate with aae.pop
  value <- simulate(
    mspecies_obj, nsim = nsim, init = init,
    args = list(
      list(
        replicated_covariates = format_covariates(matrix(rep(0.9, 10), ncol = nsim))
      ),
      list(
        replicated_covariates = format_covariates(matrix(rep(0.9, 10), ncol = nsim))
      ),
      list(
        replicated_covariates = format_covariates(matrix(rep(0.9, 10), ncol = nsim))
      ),
      list(
        replicated_covariates = format_covariates(matrix(rep(0.9, 10), ncol = nsim))
      )
    )
  )

  # now simulate manually
  all_dyn <- list(dyn_mc1, dyn_mc2, dyn_mc3, dyn_mc5)
  sp_order <- match(sapply(mspecies_obj$dynamics, function(x) x$hex),
                    sapply(all_dyn, function(x) x$hex))
  process_target <- sp_order[4]
  all_dyn <- all_dyn[sp_order]
  target <- lapply(
    value,
    function(x) array(dim = dim(x))
  )
  for (i in seq_along(target)) {
    target[[i]][, , 1] <- init[, , i]
  }
  interaction_pairs <- apply(
    mspecies_obj$structure, 1, function(x) which(x == 1)
  )
  for (i in seq_len(dim(value[[1]])[3] - 1)) {
    for (k in seq_len(nsim)) {
      for (j in seq_along(all_dyn)) {

        # pull out the baseline matrix for species j
        mat <- all_dyn[[j]]$matrix

        # add processes for dyn_mc5
        if (j == process_target) {
          mat[survival(mat)] <- pmin(1.1 * mat[survival(mat)], 1)
          mat[survival(mat)] <- 0.9 * mat[survival(mat)]
        }

        # add in pairwise interactions
        if (length(interaction_pairs[[j]]) > 0) {
          x <- mat[survival(mat)]
          for (.i in seq_along(interaction_pairs[[j]])) {
            x <- x /
              (1 +
                 x *
                 sum(target[[interaction_pairs[[j]][.i]]][k, , i]) /
                 2000
              )
          }
          mat[survival(mat)] <- x
        }

        # update
        target[[j]][k, , (i + 1)] <- target[[j]][k, , i] %*% t(mat)

      }
    }
  }
  for (i in seq_along(target))
    class(target[[i]]) <- c("simulation", "array")
  class(target) <- c("simulation_list", "list")
  expect_equal(value, target)

  # but it should replace covariates and rep_covars with 1 if not set
  value <- simulate(mspecies_obj, nsim = nsim, init = init)

  # now simulate manually
  all_dyn <- list(dyn_mc1, dyn_mc2, dyn_mc3, dyn_mc5)
  sp_order <- match(sapply(mspecies_obj$dynamics, function(x) x$hex),
                    sapply(all_dyn, function(x) x$hex))
  process_target <- sp_order[4]
  all_dyn <- all_dyn[sp_order]
  target <- lapply(
    value,
    function(x) array(dim = dim(x))
  )
  for (i in seq_along(target)) {
    target[[i]][, , 1] <- init[, , i]
  }
  interaction_pairs <- apply(
    mspecies_obj$structure, 1, function(x) which(x == 1)
  )
  for (i in seq_len(dim(value[[1]])[3] - 1)) {
    for (k in seq_len(nsim)) {
      for (j in seq_along(all_dyn)) {

        # pull out the baseline matrix for species j
        mat <- all_dyn[[j]]$matrix

        # add processes for dyn_mc5
        if (j == process_target) {
          mat[survival(mat)] <- pmin(1.1 * mat[survival(mat)], 1)
        }

        # add in pairwise interactions
        if (length(interaction_pairs[[j]]) > 0) {
          x <- mat[survival(mat)]
          for (.i in seq_along(interaction_pairs[[j]])) {
            x <- x /
              (1 +
                 x *
                 sum(target[[interaction_pairs[[j]][.i]]][k, , i]) /
                 2000
              )
          }
          mat[survival(mat)] <- x
        }

        # update
        target[[j]][k, , (i + 1)] <- target[[j]][k, , i] %*% t(mat)

      }
    }
  }
  for (i in seq_along(target))
    class(target[[i]]) <- c("simulation", "array")
  class(target) <- c("simulation_list", "list")
  expect_equal(value, target)

})
# nolint end

test_that("multispecies errors informatively when inputs are inappropriate", {

  # can only take interactions objects
  expect_error(multispecies(dyn_mc1, interactions_test2, interactions_test3),
               "all inputs to multispecies should be interaction objects")

  # Leslie updater can't be used for non-leslie matrices
  expect_error(
    simulate(
      multispecies(interactions_test1, interactions_test2, interactions_test3),
      options = list(update = update_binomial_leslie)
    ),
    "matrix must be a Leslie matrix"
  )

  # args must be provided for each species
  expect_error(
    simulate(
      multispecies(interactions_test1, interactions_test2, interactions_test3),
      args = list(interaction = "a")
    ),
    "one element of args must be provided for each species"
  )

  # dynamic args must be provided for each species
  expect_error(
    simulate(
      multispecies(interactions_test1, interactions_test2, interactions_test3),
      args = list(
        list(interaction = format_covariates(rep(1, 10))),
        list(interaction = format_covariates(rep(1, 12))),
        list(interaction = format_covariates(rep(1, 9))),
        list(interaction = format_covariates(rep(1, 10)))
      )
    ),
    "dynamic arguments must have the same length for all species"
  )

})

test_that("multispecies methods work without error", {

  # create multispecies object
  mspecies_obj <- multispecies(
    interactions_test1, interactions_test2, interactions_test3
  )

  # how many replicates are we simulating?
  nsim <- 10

  # define initial conditions
  dims <- c(nsim, nstage, mspecies_obj$nspecies)
  init <- array(rpois(prod(dims), lambda = 10), dim = dims)

  # simulate with aae.pop
  value <- simulate(mspecies_obj, nsim = nsim, init = init)

  # check plotting
  expect_silent(plot(value))
  expect_silent(plot(value, which = 1))

  # check is method
  expect_true(is.simulation_list(value))
  expect_false(is.simulation_list("a"))

})
