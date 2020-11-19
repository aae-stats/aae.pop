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
#' @param masks a logical matrix or vector (or list of these)
#'   defining cells affected by \code{funs}. See Details and
#'   \code{\link{masks}}
#' @param funs a function or list of functions with one element
#'   for each element of \code{masks}. See Details
#' @param nmask logical vector or list of vectors defining
#'   elements of the popualtion vector affected by each
#'   mask-function pair. Intended primarily for internal
#'   use when scaling up processes in
#'   \code{\link{metapopulation}}
#'
#' @details \code{density_dependence} specifies standard
#'   density dependence on vital rates, such as scramble or
#'   contest competition or allee effects.
#'
#'   \code{density_dependence_n} is an alternative
#'   parameterisation of density dependence that acts directly
#'   on population abundances.
#'
#'   Masks must be of the same dimension as the population
#'   dynamics matrix and specify cells influenced by density
#'   dependence according to \code{funs}. In the case of
#'   \code{density_dependence_n}, \code{masks} are
#'   logical vectors with one element for each class.
#'   Additional details on masks are provided
#'   in \code{\link{masks}}.
#'
#'   If using \code{density_depenence}, functions must take at
#'   least two arguments, a matrix \code{x} and a vector \code{n},
#'   which represent the population dynamics matrix and the
#'   population abundances. Functions must return a matrix with
#'   the same dimensions as \code{x}, modified to reflect the
#'   effects of current abundances (\code{n}) on
#'   vital rates.
#'
#'   In the case of \code{density_dependence_n},
#'   \code{funs} takes only one argument, the population
#'   abundances \code{n} following all other updates in a
#'   given iteration/generation. This allows rescaling of
#'   population abundances based on total abundance or
#'   through more complicated functions that depend
#'   on external arguments (e.g., mass mortality events or
#'   harvesting).
#'
#'   Additional arguments to functions are supported and can be
#'   passed to \code{\link{simulate}} with the \code{args},
#'   \code{args.dyn}, or \code{args.fn} arguments.
#'
#' @examples
#' # add
density_dependence <- function(masks, funs, nmask = NULL) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)
  force(nmask)

  if (is.list(masks)) {
    fn <- function(x, n, ...) {
      if (!is.null(nmask))
        n <- n[nmask]
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]], n, ...)
      x
    }
  } else {
    fn <- function(x, n, ...) {
      if (!is.null(nmask))
        n <- n[nmask]
      do_mask(x, masks, funs, n, ...)
    }
  }

  as_density_dependence(fn)

}

#' @rdname density_dependence
#'
#' @export
density_dependence_n <- function(masks, funs) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

  if (is.list(masks)) {
    fn <- function(pop_t, ...) {
      for (i in seq_along(masks))
        pop_t <- do_mask(pop_t, masks[[i]], funs[[i]], ...)
      pop_t
    }
  } else {
    fn <- function(pop_t, ...) {
      do_mask(pop_t, masks, funs, ...)
    }
  }

  as_density_dependence_n(fn)

}

# internal function: set density_dependence class
as_density_dependence <- function(x) {
  as_class(x, name = "density_dependence", type = "function")
}

# internal function: set density_dependence_n class
as_density_dependence_n <- function(x) {
  as_class(x, name = "density_dependence_n", type = "function")
}
