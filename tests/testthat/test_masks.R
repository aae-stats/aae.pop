context("masks")

# set up: define a population matrix
mat <- matrix(1:25, ncol = 5)

test_that("mask definitions are correct", {

  # reproduction without dims
  value <- reproduction(mat)
  target <- row(mat) == 1 & col(mat) > 1
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # reproduction with consecutive dims
  value <- reproduction(mat, dims = 4:5)
  target <- row(mat) == 1 & col(mat) > 3
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # reproduction with non-consecutive dims
  value <- reproduction(mat, dims = c(2, 5))
  target <- row(mat) == 1 & col(mat) %in% c(2, 5)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # survival without dims
  value <- survival(mat)
  target <- row(mat) == col(mat)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # survival with consecutive dims
  value <- survival(mat, dims = 3:5)
  target <- row(mat) == col(mat) & col(mat) > 2
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # survival with non-consecutive dims
  value <- survival(mat, dims = c(1, 3, 4))
  target <- row(mat) == col(mat) & col(mat) %in% c(1, 3, 4)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # transition without dims
  value <- transition(mat)
  target <- row(mat) == col(mat) + 1
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # transition with consecutive dims
  value <- transition(mat, dims = 3:4)
  target <- row(mat) == col(mat) + 1 & col(mat) > 2
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # transition with non-consecutive dims
  value <- transition(mat, dims = c(1, 3, 4))
  target <- row(mat) == col(mat) + 1 & col(mat) %in% c(1, 3, 4)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # all_cells (dims ignored)
  value <- all_cells(mat)
  target <- !is.na(mat)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # all_cells with and without dims should be identical
  target <- all_cells(mat, dims = 1)
  expect_equal(value, target)

  # all_classes without dims
  value <- all_classes(mat)
  target <- matrix(!is.na(seq_len(ncol(mat))), ncol = 1)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # all_classes with consecutive dims
  value <- all_classes(mat, dims = 4:5)
  target <- row(target) > 3
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

  # all_classes with non-consecutive dims
  value <- all_classes(mat, dims = c(2, 5))
  target <- matrix(row(value) %in% c(2, 5), ncol = 1)
  class(target) <- c("mask", "matrix", "array")
  expect_equal(value, target)

})


test_that("masks combine correctly as matrices", {

  # create a mask from two separate masks
  mask1 <- reproduction(mat, dims = c(4:5))
  mask2 <- survival(mat)
  value <- combine(mask1, mask2)
  target <- row(mat) == col(mat) |
    (row(mat) == 1 & col(mat) > 3)
  class(target) <- c("mask", "matrix", "array")
  expect_equivalent(value, target)

})

test_that("masks combine correctly as functions", {

  # create a mask function from two separate fns
  value <- combine(reproduction, survival)(mat)
  target <- row(mat) == col(mat) | row(mat) == 1
  class(target) <- c("mask", "matrix", "array")
  expect_equivalent(value, target)

})

test_that("combining masks errors informatively for inappropriate classes", {

  # errors if given non-mask/function
  expect_error(
    combine(matrix(1:10)),
    "combine is not defined for objects of class matrix"
  )

  # errors if given non-mask/function
  expect_error(
    combine(matrix(1:10), reproduction),
    "combine is not defined for objects of class matrix"
  )

  # errors if given combo including non-mask/function
  expect_error(
    combine(reproduction(mat), matrix(1:10)),
    "combine is not defined for mask objects combined with objects of class matrix"
  )

  # errors if given combo including non-mask
  expect_error(
    combine(reproduction, matrix(1:10)),
    "combine is not defined for function objects combined with objects of class matrix"
  )

})
