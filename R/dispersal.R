#' @name dispersal
#' @title Specify dispersal between populations in a metapopulation model
#' @description Specify disperal between populations, including stochasticity
#'   and density dependence in dispersal parameters
NULL

#' @rdname dispersal
#'
#' @export
#'
#' @param kernel numeric matrix specifying the probability of
#'   specific classes moving between two populations. Matrices
#'   have the same columns-move-to-rows structure as in
#'   the population dynamics matrices described in
#'   \code{\link{dynamics}}, so a non-zero value in cell (a, b)
#'   denotes a transition from class b in the source population
#'   to class a in the receiving population
#' @param stochasticity_masks a logical matrix or list of logical matrices
#'   defining cells affected by \code{stochasticity_funs}.
#'   See Details and \code{\link{masks}}
#' @param stochasticity_funs a function or list of functions with
#'   one element for each element of \code{stochasticity_masks}.
#'   See Details
#' @param density_masks a logical matrix or list of logical matrices
#'   defining cells affected by \code{density_funs}.
#'   See Details and \code{\link{masks}}
#' @param density_funs a function or list of functions with
#'   one element for each element of \code{density_masks}.
#'   See Details
#' @param proportion logical indicating whether \code{kernel}
#'   is specified in absolute probabilites or as a proportion
#'   of the source population (defaults to \code{FALSE}).
#'   If \code{TRUE}, values in \code{kernel} are calculated as
#'   a proportion of the total probability an individual
#'   exits that class at any given time step
#'
#' @details To be completed.
#'
#' @examples
#' # add
dispersal <- function(kernel,
                      stochasticity_masks = NULL,
                      stochasticity_funs = NULL,
                      density_masks = NULL,
                      density_funs = NULL,
                      proportion = FALSE) {

  # check kernel
  kernel <- check_structure(kernel)

  # define functions to calculate stochasticity for each element
  stoch_fn <- NULL
  if (!is.null(stochasticity_masks)) {
    if (is.list(stochasticity_masks)) {
      stoch_fn <- function(x) {
        for (i in seq_along(stochasticity_masks))
          x <- do_mask(x, stochasticity_masks[[i]], stochasticity_funs[[i]])
        x
      }
    } else {
      stoch_fn <- function(x) {
        do_mask(x, stochasticity_masks, stochasticity_funs)
      }
    }
  }

  # define functions to calculate density dependence for each element
  #   based on vector of abundances in the source and destination pops
  dens_fn <- NULL
  if (!is.null(density_masks)) {
    if (is.list(density_masks)) {
      dens_fn <- function(x, nsource, ndestination) {
        for (i in seq_along(density_masks))
          x <- do_mask(x, density_masks[[i]], density_funs[[i]], nsource, ndestination)
        x
      }
    } else {
      dens_fn <- function(x, nsource, ndestination) {
        do_mask(x, density_masks, density_funs, nsource, ndestination)
      }
    }
  }

  # collate output into a list
  dispersal <- list(
    kernel = kernel,
    stochasticity = stoch_fn,
    density = dens_fn,
    proportion = proportion
  )

  # and return as dispersal object
  as_dispersal(dispersal)

}

# internal function: check the structure matrix is OK
check_structure <- function(x) {

  # is structure actually a matrix?
  if (!is.matrix(x))
    stop("structure must be a matrix", call. = FALSE)

  # structure must be a square matrix
  if (nrow(x) != ncol(x))
    stop("structure must be a square matrix", call. = FALSE)

  # convert binary matrix to logical
  if (all(x %in% c(0, 1)))
    x <- x > 0

  # is structure now logical?
  if (!all(x %in% c(TRUE, FALSE)))
    stop("structure must be binary (0/1) or logical (TRUE/FALSE)", call. = FALSE)

  # remove diag if included
  if (any(diag(x)))
    diag(x) <- FALSE

  # how many populations are included?
  npop <- nrow(x)

  # how many dispersal elements exist?
  ndispersal <- sum(x[upper.tri(x)]) + sum(x[lower.tri(x)])

  # return list of key information
  list(npop = npop,
       ndispersal = ndispersal,
       structure = x)

}

# internal function: set demographic_stochasticity class
as_dispersal <- function(x) {
  as_class(x, name = "dispersal", type = "list")
}
