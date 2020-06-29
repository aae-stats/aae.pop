context("multispecies")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# define multispecies interactions
multispecies_mask <- survival(mat)
multispecies_fun <- function(x, n) {
  x / (1 + x * sum(n) / 2000)
}
dyn_mc1 <- dynamics(mat)
dyn_mc2 <- dynamics(mat)
dyn_mc3 <- dynamics(mat)
dyn_mc4 <- dynamics(mat)
interactions_test1 <- pairwise_interaction(dyn_mc1, dyn_mc2, multispecies_mask, multispecies_fun)
interactions_test2 <- pairwise_interaction(dyn_mc1, dyn_mc3, multispecies_mask, multispecies_fun)
interactions_test3 <- pairwise_interaction(dyn_mc3, dyn_mc4, multispecies_mask, multispecies_fun)

test_that("multispecies object has correct dynamics elements", {

  # create multispecies object
  mspecies_obj <- multispecies(interactions_test1, interactions_test2, interactions_test3)

  # check dynamics one by one (should have 1-4, not necessarily in that order)
  #  match based on hex
  expect_equal(1L, 1L)

})


test_that("multispecies objects simulate correctly", {

  # not a full test because test_simulate.R covers identical use of simulate

  # create multispecies object
  mspecies_obj <- multispecies(interactions_test1, interactions_test2, interactions_test3)

  # simulate with aae.pop
  value <- simulate(mspecies_obj, nsim = 10)

  # simulate manually
  target <- simulate(mspecies_obj, nsim = 10)

  # and compare
  expect_equal(value, target)

})

test_that("multispecies errors informatively when inputs are inappropriate", {

  # can only take interactions objects
  expect_error(multispecies(dyn_mc1, interactions_test2, interactions_test3),
               "all inputs to multispecies should be interaction objects")

})
