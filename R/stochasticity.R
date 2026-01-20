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
#' @returns \code{environmental_stochasticity} or
#'   \code{demographic_stochasticity} object specifying the way in which
#'   stochasticity should be included in a matrix population model;
#'   for use with \code{\link{dynamics}}
#'
#' @examples
#' # define a population matrix (columns move to rows)
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- dynamics(popmat)
#'
#' # note that there is only one trajectory now because
#' #   this simulation is deterministic.
#' #
#' # let's change that by adding some environmental stochasticity
#' envstoch <- environmental_stochasticity(
#'   masks = list(
#'     reproduction(popmat, dims = 4:5),
#'     transition(popmat)
#'   ),
#'   funs = list(
#'     \(x) rpois(n = length(x), lambda = x),
#'     \(x) runif(n = length(x), min = 0.9 * x, max = pmin(1.1 * x, 1))
#'   )
#' )
#'
#' # update the dynamics object and simulate from it
#' dyn <- update(dyn, envstoch)
#' sims <- simulate(
#'   dyn,
#'   init = c(50, 20, 10, 10, 5),
#'   nsim = 100,
#'   options = list(ntime = 50),
#' )
#'
#' # a simple way to add demographic stochasticity is to change
#' #   the "updater" that converts the population at time t
#' #   to its value at time t + 1. The default in aae.pop
#' #   uses matrix multiplication of the vital rates matrix
#' #   and current population. A simple tweak is to update
#' #   with binomial draws. Note that this also requires a
#' #   change to the "tidy_abundances" option so that population
#' #   abundances are always integer values.
#' sims <- simulate(
#'   dyn,
#'   init = c(50, 20, 10, 10, 5),
#'   nsim = 100,
#'   options = list(
#'     update = update_binomial_leslie,
#'     tidy_abundances = floor
#'   ),
#'   args = list(
#'     environmental_stochasticity = list(envstoch_function)
#'   )
#' )
environmental_stochasticity <- function(masks, funs) {
  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

  if (is.list(masks)) {
    fn <- function(x, ...) {
      for (i in seq_along(masks)) {
        x <- do_mask(x, masks[[i]], funs[[i]], ...)
      }
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
      for (i in seq_along(masks)) {
        x <- do_mask(x, masks[[i]], funs[[i]], ...)
      }
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
    x,
    name = "environmental_stochasticity", type = "function"
  )
}

# internal function: set demographic_stochasticity class
as_demographic_stochasticity <- function(x) {
  as_class(
    x,
    name = "demographic_stochasticity", type = "function"
  )
}
