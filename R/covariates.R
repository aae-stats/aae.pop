#' @name covariates
#'
#' @title Specify covariate dependence in models of population dynamics
#'
#' @description Specify relationship between a vector or matrix of covariates
#'   and vital rates.
#'
#' @export
#'
#' @param x a vector or matrix of covariates with one observation
#'   for each time step. \code{\link{dynamics}} passes \code{x} directly
#'   to \code{fun} and assumes that \code{x} is a matrix with
#'   observations in rows. Numeric vectors are converted to
#'   matrices with one column
#' @param fun a function that takes two arguments, the first
#'   being the population dynamics matrix and the second being
#'   \code{x}. Arguments can have any name but are assumed
#'   to be in this order
#'
#' @details To be completed.
#'
#' @examples
#' # add
covariates <- function(x, masks, funs) {

  # convert all x to matrix with time slices in columns
  if (!is.matrix(x))
    x <- matrix(x, ncol = 1)

  # define function to combine masks and funs
  if (is.list(masks)) {
    fn <- function(x, ...) {
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]], ...)
      x
    }
  } else {
    fn <- function(x, ...) {
      do_mask(x, masks, funs, ...)
    }
  }

  # return
  as_covariates(list(x = x, fun = fn, ntime = nrow(x)))

}

# internal function: set covariates class
as_covariates <- function(x) {
  as_class(x, name = "covariates", type = "list")
}
