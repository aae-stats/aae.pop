#' @name covariates
#'
#' @title Specify covariate dependence in models of population dynamics
#'
#' @description Specify relationship between a vector or matrix of covariates
#'   and vital rates.
#'
#' @export
#'
#' @param masks a logical matrix or vector (or list of these)
#'   defining cells affected by \code{funs}. See Details and
#'   \code{\link{masks}}
#' @param funs a function or list of functions with one element
#'   for each element of \code{masks}. See Details
#'
#' @details Masks must be of the same dimension as the population
#'   dynamics matrix and specify cells influenced by covariates
#'   according to \code{funs}. Functions must take at least
#'   one arguments, a matrix representing the population
#'   dynamics matrix. Functions must return a matrix with
#'   the same dimensions as the input, modified to reflect the
#'   effects of covariates on vital rates. Additional arguments
#'   can be passed to \code{funs} and can be specified as
#'   \code{args} in \code{\link{simulate}}.
#'
#' @examples
#' # add
covariates <- function(masks, funs) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

  # define function to combine masks and funs
  if (is.list(masks)) {
    fn <- function(mat, ...) {
      for (i in seq_along(masks))
        mat <- do_mask(mat, masks[[i]], funs[[i]], ...)
      mat
    }
  } else {
    fn <- function(mat, ...) {
      do_mask(mat, masks, funs, ...)
    }
  }

  # return
  as_covariates(fn)

}

# internal function: set covariates class
as_covariates <- function(x) {
  as_class(x, name = "covariates", type = "function")
}
