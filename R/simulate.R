#' @title Simulate single or multispecies population dynamics in R
#'
#' @name simulate
#'
#' @importFrom stats rpois
#'
#' @export
#'
#' @param obj a \code{dynamics} object created with
#'   \code{\link{define_dynamics}} or from a subsequent call to
#'   \code{\link{combine_species}} or \code{\link{expand_spatially}}
#' @param init an array of initial conditions with one row per replicate and one
#'   column per population stage. Additionally requires one slice per species if
#'   \code{obj} has been created with \code{\link{combine_species}}. Defaults
#'   to \code{NULL}, in which case initial conditions are generated randomly
#'   from a Poisson distribution
#' @param options a named \code{list} of simulation options. Currently accepted
#'   values are:
#'   - \code{ntime} the number of time steps to simulate, ignored if \code{obj}
#'       includes a \code{\link{modifier}} (default = 50)
#'   - \code{replicates} the number of replicate simulations (default = 1000)
#'   - \code{keep_slices} \code{logical} defining whether to keep intermediate
#'       population abundances or (if \code{FALSE}) to return only the final
#'       time slice
#'   - \code{tidy_abundances} a function to handle predicted abundance data
#'       that may be non-integer. Defaults to \code{identity}; suggested
#'       alternatives are \code{floor}, \code{round}, or \code{ceiling}
#'   - \code{lambda_init} lambda defining a Poisson distribution, used to
#'       generate rando initial conditions if \code{init = NULL}
#'
#' @details
#'
#' @examples
#'
#' # define a population matrix
#'
#' # define a dynamics objet
#'
#' # simulate from this
simulate <- function(obj, init = NULL, options = list()) {

  opt <- list(
    ntime = 50,
    replicates = 1000,
    keep_slices = TRUE,
    tidy_abundances = identity,
    lambda_init = 10
  )
  opt[names(options)] <- options

  if (!is.null(obj$covariates))
    opt$ntime <- obj$ntime

  pop <- initialise(obj, opt, init)

  for (i in seq_len(opt$ntime)) {

    if (obj$nspecies > 1) {
      pop[, , , i + 1] <- simulate_once_multispecies(
        obj, pop[, , , i], tidy_abundances = opt$tidy_abundances
      )
    } else {
      pop[, , i + 1] <- simulate_once(
        obj, pop[, , i], tidy_abundances = opt$tidy_abundances
      )
    }

  }

  if (!opt$keep_slices)
    pop <- pop[, , opt$ntime + 1]

  # return
  as_simulation(pop)

}

# internal function: update a single time step for one species
simulate_once <- function(obj, pop_t, tidy_abundances) {

  if (!is.null(obj$covariates)) {
    mat <- obj$matrix[[i]]
  } else {
    mat <- obj$matrix
  }

  if (!is.null(obj$environmental_stochasticity))
    mat <- obj$environmental_stochasticity(mat)

  if (!is.null(obj$density_dependence))
    mat <- obj$density_dependence(mat, pop_t)

  pop_tp1 <- tcrossprod(pop_t, mat)

  if (!is.null(obj$demographic_stochasticity))
    pop_tp1 <- t(apply(pop_tp1, 1, obj$demographic_stochasticity))

  if (!is.null(obj$density_dependence_n))
    pop_tp1 <- t(apply(pop_tp1, 1, obj$density_dependence_n))

  tidy_abundances(pop_tp1)

}

# internal function: update a single time step with interacting species
simulate_once_multispecies <- function(obj, pop_t, nspecies, tidy_abundances) {

  # calculate density effects of other species
  #   according to something save in obj (interaction matrix with fns?)

  # initialise abundance array for next time step
  pop_tp1 <- array(NA, dim = dim(pop_t))

  # loop through species, updating one-by-one
  for (i in seq_along(obj$nspecies)) {

    # apply density rescaling (don't recalculate dnesities here because will change)

    if (!is.null(obj$covariates)) {
      mat <- obj$matrix[, , i]
    } else {
      mat <- obj$matrix
    }

    if (!is.null(obj$environmental_stochasticity))
      mat <- obj$environmental_stochasticity(mat)

    if (!is.null(obj$density_dependence))
      mat <- obj$density_dependence(mat, pop_t[, , i])

    pop_tp1[, , i] <- tcrossprod(pop_t[, , i], mat)

    if (!is.null(obj$demographic_stochasticity))
      pop_tp1[, , i] <- t(apply(pop_tp1[, , i], 1, obj$demographic_stochasticity))

    if (!is.null(obj$density_dependence_n))
      pop_tp1[, , i] <- t(apply(pop_tp1[, , i], 1, obj$density_dependence_n))

  }

  tidy_abundances(pop_tp1)

}

# internal function: initialise a simulation when inits not provided
initialise <- function(obj, opt, init) {

  if (obj$nspecies > 1) {
    dims <- c(opt$replicates, obj$nclass, obj$nspecies, opt$ntime + 1)
  } else {
    dims <- c(opt$replicates, obj$nclass, opt$ntime + 1)
  }
  ndim <- length(dims)

  if (is.null(init)) {
    init <- array(
      rpois(prod(dims[-ndim]), lambda = opt$lambda_init),
      dim = dims[-ndim]
    )
  }

  pop <- array(NA, dim = dims)
  pop[seq_len(prod(dims[-ndim]))] <- init

  pop

}

# internal function: set simulation class
as_simulation <- function(x) {
  as_class(x, name = "simulation", type = "array")
}
