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
demo_fn <- function(x) rpois(length(x), lambda = x)
demostoch <- demographic_stochasticity(masks = demo_mask, funs = demo_fn)

# add environmental stochasticity
mask_list <- list(reproduction(mat, dims = 4:5), survival(mat))
fn_list <- list(
  function(x) x + rpois(length(x), lambda = 2),
  function(x) x + runif(length(x), min = -min(x), max = max(1 - x))
)
envstoch <- environmental_stochasticity(masks = mask_list, funs = fn_list)

# add rescaling method (alternative form of density dependence)
rescale_mask <- all_stages(mat)
rescale_fn <- function(x) 200 * (x / sum(x))
resc <- density_dependence_n(rescale_mask, rescale_fn)

# manual simulation function
sim_manual <- function(mat, ntime, nsim, init) {
  nstage <- ncol(mat)
  out <- array(NA, dim = c(nsim, nstage, ntime + 1))
  out[, , 1] <- init
  for (i in seq_len(ntime))
    out[, , (i + 1)] <- floor(out[, , i] %*% t(mat))
  out
}

test_that("simulate returns correct abundances with different processes", {

  # simulate trajectories with no uncertainty
  nsim <- 10
  ntime <- 50
  dyn <- dynamics(mat)
  init_set <- matrix(rpois(nstage * nsim, lambda = 20), ncol = nstage)
  value <- simulate(
    dyn,
    nsim = nsim,
    init = init_set,
    options = list(ntime = ntime, tidy_abundances = floor)
  )
  target <- sim_manual(
    mat, ntime = ntime, nsim = nsim, init = init_set
  )
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
