#' @name dispersal
#' @title Specify dispersal between populations in a metapopulation model
#' @description Specify dispersal between populations, including stochasticity
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
#'   is specified in absolute probabilities or as a proportion
#'   of the source population (defaults to \code{FALSE}).
#'   If \code{TRUE}, values in \code{kernel} are calculated as
#'   a proportion of the total probability an individual
#'   exits that class at any given time step
#'
#' @details To be completed.
#'
#' @examples
#' # define some populations, all with identical vital rates
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- lapply(
#'   1:3,
#'   function(i) dynamics(popmat)
#' )
#'
#' # define metapopulation structure with populations
#' #   1 and 3 dispersing into population 2
#' pop_structure <- matrix(0, nrow = 3, ncol = 3)
#' pop_structure[1, 2] <- 1
#' pop_structure[3, 2] <- 1
#'
#' # define dispersal between populations
#' dispersal_matrix <- matrix(0, nrow = nclass, ncol = nclass)
#' dispersal_matrix[survival(dispersal_matrix, dims = 20:25)] <- 0.2
#' pop_dispersal1 <- dispersal(dispersal_matrix, proportion = TRUE)
#' pop_dispersal2 <- dispersal(dispersal_matrix, proportion = FALSE)
#' pop_dispersal <- list(pop_dispersal1, pop_dispersal2)
#'
#' # create metapopulation object
#' metapop <- metapopulation(pop_structure, dyn, pop_dispersal)
#'
#' # simulate without covariates
#' sims <- simulate(metapop, nsim = 10)
#'
#' # and plot the simulated trajectories
#' plot(sims)
dispersal <- function(kernel,
                      stochasticity_masks = NULL,
                      stochasticity_funs = NULL,
                      density_masks = NULL,
                      density_funs = NULL,
                      proportion = FALSE) {
  # check kernel
  if (nrow(kernel) != ncol(kernel)) {
    stop("kernel must be a square matrix", call. = FALSE)
  }

  # force evaluation to avoid NULL functions down the line
  force(stochasticity_masks)
  force(stochasticity_funs)
  force(density_masks)
  force(density_funs)

  # define functions to calculate stochasticity for each element
  stoch_fn <- NULL
  if (!is.null(stochasticity_masks)) {
    if (is.list(stochasticity_masks)) {
      stoch_fn <- function(x, ...) {
        for (i in seq_along(stochasticity_masks)) {
          x <- do_mask(
            x, stochasticity_masks[[i]], stochasticity_funs[[i]], ...
          )
        }
        x
      }
    } else {
      stoch_fn <- function(x, ...) {
        do_mask(x, stochasticity_masks, stochasticity_funs, ...)
      }
    }
  }

  # define functions to calculate density dependence for each element
  #   based on vector of abundances in the source and destination pops
  dens_fn <- NULL
  if (!is.null(density_masks)) {
    if (is.list(density_masks)) {
      dens_fn <- function(x, n, ...) {
        for (i in seq_along(density_masks)) {
          x <- do_mask(
            x,
            density_masks[[i]],
            density_funs[[i]],
            n,
            ...
          )
        }
        x
      }
    } else {
      dens_fn <- function(x, n, ...) {
        do_mask(
          x, density_masks, density_funs, n, ...
        )
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

# internal function: set demographic_stochasticity class
as_dispersal <- function(x) {
  as_class(x, name = "dispersal", type = "list")
}
