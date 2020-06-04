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

# internal function: set demographic_stochasticity class
as_dispersal <- function(x) {
  as_class(x, name = "dispersal", type = "list")
}
