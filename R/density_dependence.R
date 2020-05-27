#' @name density_dependence
#' @title Specify density dependence in models of population dynamics
#' @description Specify density dependence in vital rates
#'   (\code{density_dependence}) and in total abundances
#'   (\code{density_dependence_n}).
NULL

#' @rdname density_dependence
#'
#' @export
#'
#' @param masks a logical matrix or list of logical matrices
#'   defining cells affected by \code{funs}. See Details and
#'   \code{\link{masks}}
#' @param funs a function or list of functions with one element
#'   for each element of \code{masks}. See Details
#'
#' @details Masks must be of the same dimension as the population
#'   dynamics matrix and specify cells influenced by density
#'   dependence according to \code{funs}. Functions must take two
#'   arguments, a matrix \code{x} and a vector \code{n}, which
#'   represent the population dynamics matrix and the population
#'   abundances. Functions must return a matrix with
#'   the same dimensions as \code{x}, modified to reflect the
#'   effects of current abundances by class (\code{n}) on
#'   vital rates
#'
#' @examples
#' # add
density_dependence <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(x, n) {
      for (i in seq_along(masks))
        x[masks[[i]]] <- funs[[i]](x[masks[[i]]], n)
      x
    }
  } else {
    fn <- function(x, n) {
        x[masks] <- funs(x[masks], n)
        x
    }
  }

  as_density_dependence(fn)

}

#' @rdname density_dependence
#'
#' @export
#'
#' @details something
#'
#' @examples
#' # add
density_dependence_n <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(pop_t) {
      for (i in seq_along(masks))
        pop_t[masks[[i]]] <- funs[[i]](pop_t[masks[[i]]])
      pop_t
    }
  } else {
    fn <- function(pop_t) {
      pop_t[masks] <- funs(pop_t[masks])
      pop_t
    }
  }

  as_density_dependence_n(fn)

}

#' @rdname density_dependence
#'
#' @export
#'
#' @param k carrying capacity used to define models of
#'   density dependence. See details for details of
#'   currently implemented models and their parameters.
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence. Currently implemented models
#'   are the Ricker model and Beverton-Holt model, both with
#'   a single parameter K.
#'
#' @examples
#' # add
beverton_holt <- function(k) {

  function(x, n) {
    x / (1 + x * sum(n) / k)
  }

}

#' @rdname density_dependence
#'
#' @export
#'
#' @examples
#' # add
ricker <- function(k) {
  NULL
}

# internal function: set density_dependence class
as_density_dependence <- function(x) {
  as_class(x, name = "density_dependence", type = "function")
}

# internal function: set density_dependence_n class
as_density_dependence_n <- function(x) {
  as_class(x, name = "density_dependence_n", type = "function")
}
