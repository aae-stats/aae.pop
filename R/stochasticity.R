#' @name stochasticity
#' @title Specify environmental and demographic stochasticity in models of
#'   population dynamics
#' @description Specify environmental stochasticity (random variation
#'   in vital rates) and demographic stochasticity (random variation
#'   in population outcomes).
NULL

#' @rdname stochasticity
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
#'   dynamics matrix (in the case of environmental stochasticity)
#'   or have one element for each class (in the case of demographic
#'   stochasticity). Masks specify cells influenced by stochasticity
#'   according to \code{funs}. Additional details on masks are provided
#'   in \code{\link{masks}}.
#'
#'   Functions must have at least one argument, a population
#'   dynamics matrix for environmental
#'   stochasticity or a vector of population abundances for
#'   demographic stochasticity. Functions must return an
#'   output of the same dimensions as the input,
#'   modified to reflect the effects of stochasticity on
#'   vital rates or population abundances.
#'
#'   Additional arguments to functions are supported and can be
#'   passed to \code{\link{simulate}} with the \code{args},
#'   \code{args.dyn}, or \code{args.fn} arguments.
#'
#' @examples
#' # add
environmental_stochasticity <- function(masks, funs) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

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

  as_environmental_stochasticity(fn)

}

#' @rdname stochasticity
#'
#' @export
demographic_stochasticity <- function(masks, funs) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

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

  as_demographic_stochasticity(fn)

}

# internal function: set environmental_stochasticity class
as_environmental_stochasticity <- function(x) {
  as_class(
    x, name = "environmental_stochasticity", type = "function"
  )
}

# internal function: set demographic_stochasticity class
as_demographic_stochasticity <- function(x) {
  as_class(
    x, name = "demographic_stochasticity", type = "function"
  )
}
