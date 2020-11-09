context("metapopulation")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# define metapopulation structure (5x5 and 3x3 versions)
structure_obj5 <- matrix(0, nrow = 5, ncol = 5)
structure_obj3 <- matrix(0, nrow = 3, ncol = 3)
structure_obj5[1, 2] <- 1
structure_obj5[4, 2] <- 1
structure_obj3[1, 2] <- 1
structure_obj3[3, 2] <- 1

# define a basic dispersal matrix (20% of stages 4 and 5 moving)
dispersal_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
dispersal_mat[survival(dispersal_mat, dims = 4:5)] <- 0.2
dispersal_obj <- dispersal(dispersal_mat, proportion = TRUE)
dispersal_list <- list(dispersal_obj, dispersal_obj)

# define dynamics objects for all populations
dynamics_obj <- dynamics(mat)
dynamics_list <- list(
  dynamics_obj, dynamics_obj, dynamics_obj, dynamics_obj, dynamics_obj
)

test_that("metapopulation returns correct output matrix", {

  # define a metapopulation
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list)

  # check that diagonal band is correct
  for (i in seq_len(metapop_test$npopulation)) {
    idx <- 1:5 + (i - 1) * 5
    mat_compare <- mat
    if (i == 2)
      mat_compare[4:5, 4:5] <- 0.6 * mat_compare[4:5, 4:5]
    expect_equal(metapop_test$matrix[idx, idx], mat_compare)
  }

  # and check dispersal bits (2 to 1)
  idx <- 1:5
  idy <- 6:10
  mat_compare <- matrix(0, nrow = 5, ncol = 5)
  mat_compare[4, 4] <- 0.2 * sum(mat[4:5, 4])
  mat_compare[5, 5] <- 0.2 * mat[5, 5]
  expect_equal(metapop_test$matrix[idx, idy], mat_compare)

  # now check 2 to 4
  idx <- 16:20
  idy <- 6:10
  mat_compare <- matrix(0, nrow = 5, ncol = 5)
  mat_compare[4, 4] <- 0.2 * sum(mat[4:5, 4])
  mat_compare[5, 5] <- 0.2 * mat[5, 5]
  expect_equal(metapop_test$matrix[idx, idy], mat_compare)

  # and everything else should be zero
  mat_remain <- metapop_test$matrix
  for (i in seq_len(metapop_test$npopulation)) {
    idx <- 1:5 + (i - 1) * 5
    mat_remain[idx, idx] <- 0
  }
  mat_remain[1:5, 6:10] <- 0
  mat_remain[16:20, 6:10] <- 0
  expect_equal(mat_remain, matrix(0, nrow = 25, ncol = 25))

  # repeat main check with 3-pop metapopulation
  metapop_test <- metapopulation(
    structure_obj3, dynamics_list[1:3], dispersal_list
  )

  # check that diagonal band is correct
  for (i in seq_len(metapop_test$npopulation)) {
    idx <- 1:5 + (i - 1) * 5
    mat_compare <- mat
    if (i == 2)
      mat_compare[4:5, 4:5] <- 0.6 * mat_compare[4:5, 4:5]
    expect_equal(metapop_test$matrix[idx, idx], mat_compare)
  }

  # and check dispersal bits (2 to 1)
  idx <- 1:5
  idy <- 6:10
  mat_compare <- matrix(0, nrow = 5, ncol = 5)
  mat_compare[4, 4] <- 0.2 * sum(mat[4:5, 4])
  mat_compare[5, 5] <- 0.2 * mat[5, 5]
  expect_equal(metapop_test$matrix[idx, idy], mat_compare)

  # now check 2 to 3
  idx <- 11:15
  idy <- 6:10
  mat_compare <- matrix(0, nrow = 5, ncol = 5)
  mat_compare[4, 4] <- 0.2 * sum(mat[4:5, 4])
  mat_compare[5, 5] <- 0.2 * mat[5, 5]
  expect_equal(metapop_test$matrix[idx, idy], mat_compare)

  # and everything else should be zero
  mat_remain <- metapop_test$matrix
  for (i in seq_len(metapop_test$npopulation)) {
    idx <- 1:5 + (i - 1) * 5
    mat_remain[idx, idx] <- 0
  }
  mat_remain[1:5, 6:10] <- 0
  mat_remain[11:15, 6:10] <- 0
  expect_equal(mat_remain, matrix(0, nrow = 15, ncol = 15))

})


test_that("metapopulation objects simulate correctly", {

  # not a full test because test_simulate.R covers identical use of simulate

  # set initial conditions
  init_set <- matrix(rpois(250, 10), ncol = 25)

  # simulate with aae.pop functions
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list)
  value <- simulate(
    metapop_test, nsim = 10, init = init_set, options = list(ntime = 10)
  )

  # simulate manually but with metapopulation$matrix
  target <- array(NA, dim = c(10, 25, 11))
  target[, , 1] <- init_set
  for (i in seq_len(10))
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_test$matrix)
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

  # simulate manually with completely re-defined matrix
  metapop_matrix <- matrix(0, nrow = 25, ncol = 25)
  for (i in seq_len(5)) {
    idx <- 1:5 + (i - 1) * 5
    metapop_matrix[idx, idx] <- mat
  }
  metapop_matrix[9:10, 9:10] <- 0.6 * metapop_matrix[9:10, 9:10]
  dispersal_elements <- matrix(0, nrow = 5, ncol = 5)
  dispersal_elements[4, 4] <- 0.2 * sum(mat[4:5, 4])
  dispersal_elements[5, 5] <- 0.2 * mat[5, 5]
  metapop_matrix[1:5, 6:10] <- dispersal_elements
  metapop_matrix[16:20, 6:10] <- dispersal_elements
  target <- array(NA, dim = c(10, 25, 11))
  target[, , 1] <- init_set
  for (i in seq_len(10))
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_matrix)
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

})

test_that("metapopulation errors informatively when
           structure or dispersal are inappropriate", {

  # dispersal must have one element for each dispersal marked in structure
  expect_error(
    metapopulation(
      structure_obj5, dynamics_list, dispersal_list[1:4]
    ),
    "dispersal must have one element for each"
  )

  # dispersal must have one element for each dispersal marked in structure
  expect_error(
    metapopulation(structure_obj5, dynamics_list, dispersal_list[1]),
    "dispersal must have one element for each"
  )

  # dims of structure must match number of dynamics objects
  expect_error(
    metapopulation(structure_obj3, dynamics_list[1:5], dispersal_list),
    "dynamics must be a single dynamics object or a list of dynamics"
  )

  # structure isn't a matrix
  expect_error(
    metapopulation(
      list(structure_obj3), dynamics_list[1:3], dispersal_list
    ),
    "structure must be a matrix"
  )

  # structure isn't a square matrix
  expect_error(
    metapopulation(
      structure_obj3[1:2, ], dynamics_list[1:3], dispersal_list
    ),
    "structure must be a square matrix"
  )

})
