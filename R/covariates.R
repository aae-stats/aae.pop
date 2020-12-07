#' @name covariates
#' @title Specify covariate dependence in models of population dynamics
#' @description Specify relationship between a vector or matrix of covariates
#'   and vital rates.
NULL

#' @rdname covariates
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
#'   according to \code{funs}. Additional details on masks are
#'   provided in \code{\link{masks}}.
#'
#'   Functions must take at least one argument, a vector or matrix
#'   representing the masked elements of the population dynamics
#'   matrix. Incorporating covariate values requires a second
#'   argument. Functions must return a vector or matrix with the
#'   same dimensions as the input, modified to reflect the
#'   effects of covariates on vital rates.
#'
#'   Additional arguments to functions are supported and can be
#'   passed to \code{\link{simulate}} with the \code{args},
#'   \code{args.dyn}, or \code{args.fn} arguments.
#'
#'   \code{format_covariates} is a helper function
#'   that takes covariates and auxiliary values as inputs and
#'   returns a correctly formatted list that can be passed
#'   as \code{args} to \code{simulate}.
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

#' @rdname covariates
#'
#' @export
#'
#' @param x a vector, matrix, or data.frame of time-varying
#'   covariate values with one element or row per time step
#' @param aux additional, static arguments to be passed to
#'   a covariates function
#' @param names optional vector of names for each covariate
#'   included in \code{x}
#'
format_covariates <- function(x, aux = NULL, names = NULL) {

  # is x a vector or matrix/data.frame?
  if (is.null(dim(x)))
    x <- matrix(x, ncol = 1)

  # add names if required
  if (!is.null(names)) {
    if (length(names) != ncol(x)) {
      stop("names must contain one value for each column of x",
           call. = FALSE)
    }
    colnames(x) <- names
  }

  # turn x into a list with one element per row/timestep
  x <- list(lapply(x, function(.x) .x[i, ]))

  # do we need to add auxiliary variables?
  if (!is.null(aux)) {
    if (!is.list(aux))
      aux <- list(aux)
    x <- c(x, aux)
  }

  # return
  x

}

# internal function: set covariates class
as_covariates <- function(x) {
  as_class(x, name = "covariates", type = "function")
}
