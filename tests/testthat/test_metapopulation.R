context("metapopulation")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- c(22, 17)
mat[survival(mat)] <- c(0.4, 0.2, 0.6, 0.6, 0.2)
mat[transition(mat)] <- c(0.56, 0.65, 0.42, 0.57)

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

# define stochasticity on dispersal
stoch_mask <- dispersal_mat > 0
stoch_fun <- function(x) {
  pmin(1.1 * x, 1)
}
dispersal_stoch <- dispersal(
  kernel = dispersal_mat,
  stochasticity_masks = stoch_mask,
  stochasticity_funs = stoch_fun,
  proportion = TRUE
)
dispersal_stoch2 <- dispersal(
  kernel = dispersal_mat,
  stochasticity_masks = list(stoch_mask, stoch_mask),
  stochasticity_funs = list(stoch_fun, stoch_fun),
  proportion = TRUE
)
dispersal_list2 <- list(dispersal_stoch, dispersal_stoch)
dispersal_list3 <- list(dispersal_stoch2, dispersal_stoch2)

# define some density terms
dd_fun1 <- function(x, n) {
  # dispersing from population 2 to 1
  x * exp(1 - sum(n[1:5]) / 20) / exp(1)
}
dd_fun2 <- function(x, n) {
  # dispersing from population 2 to 4
  x * exp(1 - sum(n[16:20]) / 20) / exp(1)
}
dispersal_dd1 <- dispersal(
  kernel = dispersal_mat,
  density_masks = stoch_mask,
  density_funs = dd_fun1,
  proportion = TRUE
)
dispersal_dd2 <- dispersal(
  kernel = dispersal_mat,
  density_masks = stoch_mask,
  density_funs = dd_fun2,
  proportion = TRUE
)
dispersal_dd3 <- dispersal(
  kernel = dispersal_mat,
  density_masks = list(stoch_mask, stoch_mask),
  density_funs = list(dd_fun1, dd_fun1),
  proportion = TRUE
)
dispersal_dd4 <- dispersal(
  kernel = dispersal_mat,
  density_masks = list(stoch_mask, stoch_mask),
  density_funs = list(dd_fun2, dd_fun2),
  proportion = TRUE
)

# combine the dispersal objects for both populations
dispersal_list4 <- list(dispersal_dd1, dispersal_dd2)
dispersal_list5 <- list(dispersal_dd3, dispersal_dd4)

# define dynamics objects for all populations
dynamics_obj <- dynamics(mat)
dynamics_list <- list(
  dynamics_obj, dynamics_obj, dynamics_obj, dynamics_obj, dynamics_obj
)

# create a mega dynamics object with all other processes in pop1
dynamics_obj2 <- update(
  dynamics_obj,
  density_dependence(
    masks = reproduction(mat, dims = 5),
    funs = \(x, n, ...) x - 0.1
  ),
  covariates(
    masks = transition(mat, dims = 3),
    funs = \(x, y, ...) x + y
  ),
  replicated_covariates(
    masks = transition(mat, dims = 4),
    funs = \(x, y, ...) x + y
  ),
  environmental_stochasticity(
    masks = reproduction(mat, dims = 4),
    funs = \(x, n, ...) x + 0.5
  ),
  demographic_stochasticity(
    masks = all_classes(mat, dims = 5),
    funs = \(x, ...) x - 1
  ),
  density_dependence_n(
    masks = all_classes(mat, dims = 4),
    funs = \(x, ...) x + 1
  ),
  add_remove_pre(
    masks = all_classes(mat, dims = 3),
    funs = \(x, ...) x + 1
  ),
  add_remove_post(
    masks = all_classes(mat, dims = 5),
    funs = \(x, ...) x + 1
  )
)
dynamics_list2 <- list(
  dynamics_obj2, dynamics_obj, dynamics_obj, dynamics_obj, dynamics_obj
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


test_that("metapopulation objects simulate correctly with processes", {

  # set initial conditions
  init_set <- matrix(rpois(25, 1), ncol = 25)

  # simulate with aae.pop functions
  metapop_test <- metapopulation(structure_obj5, dynamics_list2, dispersal_list)
  value <- simulate(
    metapop_test, nsim = 1, init = init_set, options = list(ntime = 10),
    args = list(
      covariates = list(y = 0.01),
      replicated_covariates = format_covariates(matrix(rep(0.02, 10), ncol = 1))
    )
  )

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
  metapop_matrix[1, 5] <- metapop_matrix[1, 5] - 0.1
  metapop_matrix[1, 4] <- metapop_matrix[1, 4] + 0.5
  metapop_matrix[4, 3] <- metapop_matrix[4, 3] + 0.01
  metapop_matrix[5, 4] <- metapop_matrix[5, 4] + 0.02
  target <- array(NA, dim = c(1, 25, 11))
  target[, , 1] <- init_set
  for (i in seq_len(10)) {
    tmp <- target
    tmp[1, 3, i] <- tmp[1, 3, i] + 1
    target[1, , i + 1] <- tcrossprod(tmp[1, , i], metapop_matrix)
    target[1, 5, i + 1] <- target[1, 5, i + 1] - 1
    target[1, 4, i + 1] <- target[1, 4, i + 1] + 1
    target[1, 5, i + 1] <- target[1, 5, i + 1] + 1
  }
  target[, , 1] <- init_set
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

})

test_that("metapopulation stochasticity simulates correctly", {

  # set initial conditions
  init_set <- matrix(rpois(250, 10), ncol = 25)

  # simulate with aae.pop functions
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list2)
  value <- simulate(
    metapop_test, nsim = 10, init = init_set, options = list(ntime = 10)
  )

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
  metapop_matrix[1:5, 6:10] <- pmin(1.1 * dispersal_elements, 1)
  metapop_matrix[16:20, 6:10] <- pmin(1.1 * dispersal_elements, 1)
  target <- array(NA, dim = c(10, 25, 11))
  target[, , 1] <- init_set
  for (i in seq_len(10))
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_matrix)
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

  # repeat with a list of disp stoch terms
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list3)
  value <- simulate(
    metapop_test, nsim = 10, init = init_set, options = list(ntime = 10)
  )

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
  metapop_matrix[1:5, 6:10] <- pmin(1.1 * 1.1 * dispersal_elements, 1)
  metapop_matrix[16:20, 6:10] <- pmin(1.1 * 1.1 * dispersal_elements, 1)
  target <- array(NA, dim = c(10, 25, 11))
  target[, , 1] <- init_set
  for (i in seq_len(10))
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_matrix)
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

})

test_that("metapopulation density dependence simulates correctly", {

  # set initial conditions
  init_set <- matrix(rpois(25, 10), ncol = 25)

  # simulate with aae.pop functions
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list4)
  value <- simulate(
    metapop_test, nsim = 1, init = init_set, options = list(ntime = 10)
  )

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
  target <- array(NA, dim = c(1, 25, 11))
  target[1, , 1] <- init_set
  for (i in seq_len(10)) {
    scale_fac1 <- exp(1 - sum(target[1, 1:5, i]) / 20) / exp(1)
    scale_fac2 <- exp(1 - sum(target[1, 16:20, i]) / 20) / exp(1)
    metapop_matrix[1:5, 6:10] <- dispersal_elements * scale_fac1
    metapop_matrix[16:20, 6:10] <- dispersal_elements * scale_fac2
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_matrix)
  }
  class(target) <- c("simulation", "array")

  # and compare
  expect_equal(target, value)

  # repeat with a list of dd stoch terms
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list5)
  value <- simulate(
    metapop_test, nsim = 1, init = init_set, options = list(ntime = 10)
  )
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
  target <- array(NA, dim = c(1, 25, 11))
  target[1, , 1] <- init_set
  for (i in seq_len(10)) {
    scale_fac1 <- exp(1 - sum(target[1, 1:5, i]) / 20) / exp(1)
    scale_fac2 <- exp(1 - sum(target[1, 16:20, i]) / 20) / exp(1)
    metapop_matrix[1:5, 6:10] <- dispersal_elements * scale_fac1 * scale_fac1
    metapop_matrix[16:20, 6:10] <- dispersal_elements * scale_fac2 * scale_fac2
    target[, , i + 1] <- tcrossprod(target[, , i], metapop_matrix)
  }
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

             # dispersal must exist
             expect_error(
               metapopulation(structure_obj5, dynamics_list),
               "dispersal must be provided"
             )

             # dims of populations must all match
             dynamics_tmp <- dynamics_list
             dynamics_tmp[[5]] <- dynamics(mat[1:4, 1:4])
             expect_error(
               metapopulation(structure_obj5, dynamics_tmp, dispersal_list),
               "all populations in dynamics must have the same"
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

             # structure isn't a binary matrix
             expect_error(
               metapopulation(
                 structure_obj3 * 0.5, dynamics_list[1:3], dispersal_list
               ),
               "structure must be binary"
             )

             # structure diagnonal should be overwritten
             str_test <- structure_obj3
             str_test[1, 1] <- TRUE
             expect_message(
               metapopulation(
                 str_test, dynamics_list[1:3], dispersal_list
               ),
               "exceeds 1 for classes"
             )

             # should work with a single model though
             expect_message(
               metapopulation(
                 structure_obj5, dynamics_list[[1]], dispersal_list
               ),
               "exceeds 1 for classes"
             )

           })


test_that("metapopulation S3 methods work correctly", {

  # create a metapopulation
  metapop_test <- metapopulation(structure_obj5, dynamics_list, dispersal_list4)

  # print method
  expect_output(print(metapop_test), "Metapopulation dynamics object")

  # is method
  expect_true(is.metapopulation(metapop_test))
  expect_false(is.metapopulation("a"))

})
