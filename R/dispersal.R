#' @name dispersal
#' @title Specify dispersal between populations in a metapopulation model
#' @description Specify disperal between populations, including stochasticity
#'   and density dependence in dispersal parameters
NULL

#' @rdname dispersal
#'
#' @export
#'
#' @param kernel djkfd
#' @param stochasticity_masks dfds
#' @param stochasticity_funs dfds
#' @param density_masks dfds
#' @param density_funs dfds
#'
#' @details something
#'
#' @examples
#' # add
dispersal <- function(kernel,
                      stochasticity_masks = NULL,
                      stochasticity_funs = NULL,
                      density_masks = NULL,
                      density_funs = NULL) {

  # check kernel
  if (!is.matrix(kernel))
    stop("kernel must be a two dimensional matrix", call. = FALSE)

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
    density = dens_fn
  )

  # and return as dispersal object
  as_dispersal(dispersal)

}

#' @rdname dispersal
#'
#' @export
#'
#' @param kernel kernel defined as a mask of proportions
#' @param source population dynamics matrix for source population
#'
#' @details \code{rescale_dispersal} is a helper function to assist
#'   with definitions of kernels and populations when dealing with
#'   probabilities of dispersal. This function takes a kernel
#'   of proportions and a population dynamics matrix for a source
#'   population and automatically rescales both so that a fixed
#'   proportion of the population leaves and/or remains in the
#'   population at any time step. This function aims to avoid issues
#'   where dispersal kernels are inappopriately defined, resulting
#'   in unexpected creation or destruction of individuals at each
#'   time step.
rescale_dispersal <- function(kernel, source) {

  # check dims match
  if (!all.equal(dim(kernel), dim(source)))
    stop("kernel and source must have identical dimensions", call. = FALSE)

  # calculate total probabilities in source, excluding fecundity
  idx <- options()$aae.pop_reproduction_mask(source)
  reprod <- source[idx]
  source[idx] <- 0
  leave <- apply(source, 2, sum)

  # work out proportion staying in source
  remain <- 1 - apply(kernel, 2, sum)

  # work out probability of leaving
  kernel <- sweep(kernel, 2, leave, "*")

  # and of remaining
  source <- sweep(source, 2, remain, "*")

  # add fecundity back in
  source[idx] <- reprod

  # and return both
  list(kernel = kernel,
       source = source)

}

# internal function: set demographic_stochasticity class
as_dispersal <- function(x) {
  as_class(x, name = "dispersal", type = "list")
}
