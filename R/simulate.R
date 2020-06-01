#' @title Simulate single or multispecies population dynamics in R
#'
#' @importFrom stats rpois rbinom runif simulate
#'
#' @export
#'
#' @param object a \code{dynamics} object created with
#'   \code{\link{dynamics}} or from a subsequent call to
#'   \code{\link{multispecies}} or \code{\link{metapopulation}}
#' @param nsim the number of replicate simulations (default = 1)
#' @param seed optional seed used prior to initialisation and simulation to
#'   give reproducible results
#' @param init an array of initial conditions with one row per replicate and one
#'   column per population stage. Additionally requires one slice per species if
#'   \code{obj} has been created with \code{\link{multispecies}}. Defaults
#'   to \code{NULL}, in which case initial conditions are generated randomly
#'   according to \code{options()$aae.pop_initialisation}
#' @param options a named \code{list} of simulation options. Currently accepted
#'   values are:
#'   - \code{ntime} the number of time steps to simulate, ignored if \code{obj}
#'       includes a \code{\link{covariates}} (default = 50)
#'   - \code{keep_slices} \code{logical} defining whether to keep intermediate
#'       population abundances or (if \code{FALSE}) to return only the final
#'       time slice
#'   - \code{tidy_abundances} a function to handle predicted abundance data
#'       that may be non-integer. Defaults to \code{identity}; suggested
#'       alternatives are \code{floor}, \code{round}, or \code{ceiling}
#'   - \code{initialise_args} a list of arguments passed to the function
#'       used to initialise abundance trajectories. Only used if
#'       \code{init = NULL}. Defaults to \code{options()$aae.pop_lambda},
#'       which specifies lambda for Poisson random draws. The default
#'       initialisation function is defined by
#'       \code{options()$aae.pop_initialisation}.
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
# nolint start
simulate.dynamics <- function(object,
                              nsim = 1,
                              seed = NULL,
                              init = NULL,
                              options = list()) {
  # nolint end

  # set default options for simulation
  opt <- list(
    ntime = options()$aae.pop_ntime,
    keep_slices = options()$aae.pop_keep_slices,
    tidy_abundances = options()$aae.pop_tidy_abundances,
    initialise_args = list(options()$aae.pop_lambda)
  )
  opt[names(options)] <- options

  # add nsim into options
  opt$replicates <- nsim

  # use the number of covariate values instead of fixed ntime if
  #   covariates are provided
  if (!is.null(object$covariates))
    opt$ntime <- object$ntime

  # if seed is provided, use it but reset random seed afterwards
  if (!is.null(seed)) {

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)
    }

    r_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", r_seed, envir = .GlobalEnv))
    set.seed(seed)

  }

  # initalise the population with init if provided, following
  #   options()$aae.pop_initialisation otherwise
  pop <- initialise(object, opt, init)

  # loop through timesteps, updating population at each timestep
  for (i in seq_len(opt$ntime)) {

    # split based on multispecies or single species
    if (object$nspecies > 1) {
      pop[, , , i + 1] <- simulate_once_multispecies(
        iter = i,
        object,
        pop[, , , i],
        opt = opt
      )
    } else {
      pop[, , i + 1] <- simulate_once(
        iter = i,
        object,
        pop[, , i],
        opt = opt
      )
    }

  }

  # do we want to keep intermediate abundances or just the final step?
  if (!opt$keep_slices)
    pop <- pop[, , opt$ntime + 1]

  # return
  as_simulation(pop)

}

# internal function: update a single time step for one species
simulate_once <- function(iter, obj, pop_t, opt) {

  # matrix will be a list if expanded over covariates
  if (!is.null(obj$covariates)) {
    mat <- obj$matrix[[iter]]
  } else {
    mat <- obj$matrix
  }

  # keep pop_t as a matrix if replicates == 1
  if (opt$replicates == 1)
    pop_t <- matrix(pop_t, nrow = 1)

  # draw stochastic matrix values if env stoch included
  if (!is.null(obj$environmental_stochasticity))
    mat <- obj$environmental_stochasticity(mat)

  # tweak matrix to account for density effects on vital rates,
  #   setting a flag to change update step accordingly
  is_expanded <- FALSE
  if (!is.null(obj$density_dependence)) {
    mat_list <- lapply(
      seq_len(opt$replicates),
      function(i) obj$density_dependence(mat, pop_t[i, ])
    )
    is_expanded <- TRUE
  }

  # single-step update of abundances
  if (is_expanded) {
    pop_tp1 <- t(mapply(
      options()$aae.pop_update,
      lapply(seq_len(opt$replicates), function(i) pop_t[i, ]),
      mat_list
    ))
  } else {
    pop_tp1 <- options()$aae.pop_update(pop_t, mat)
  }

  # tweak abundances to add stochastic variation in demographic
  #   outcomes
  if (!is.null(obj$demographic_stochasticity))
    pop_tp1 <- t(apply(pop_tp1, 1, obj$demographic_stochasticity))

  # final hit to abundances if they are rescaled based on biomass
  #   constraints or similar
  if (!is.null(obj$density_dependence_n))
    pop_tp1 <- t(apply(pop_tp1, 1, obj$density_dependence_n))

  # return tidied abundances (e.g. rounded or floored values)
  opt$tidy_abundances(pop_tp1)

}

# internal function: update a single time step with interacting species
simulate_once_multispecies <- function(iter,
                                       obj,
                                       pop_t,
                                       nspecies,
                                       opt) {

  # calculate density effects of other species
  #   according to something save in obj (interaction matrix with fns?)

  # initialise abundance array for next time step
  pop_tp1 <- array(NA, dim = dim(pop_t))

  # loop through species, updating one-by-one
  for (j in seq_along(obj$nspecies)) {

    # apply density rescaling
    #  (don't recalculate dnesities here because will change
    #   with each species update)

    # update as above, accounting for DIMS of pop_t carefully


  }

  opt$tidy_abundances(pop_tp1)

}

# internal function: update abundances for one time step
update_crossprod <- function(pop, mat) {
  tcrossprod(pop, mat)
}

# internal function: update abundances with a direct RNG draw
#   that combines update with demographic stochasticity
update_binomial <- function(pop, mat) {

  stop("update_binomial is not implemented", call. = FALSE)

  # check that counts are round values, otherwise can't work with size of rbinom
  # perhaps just check options()$tidy_abundances and error/warn if needed?

  # counts needs to be numbers in classes 1:(nstage-1),
  #  with nstage count added to last one
  #  not quite -- needs to be counts in classes aligned with rows??

  # surv_vec needs to be survival summed over all ways to get into a class
  surv_vec <- NULL
  fec_vec <- NULL
  counts <- NULL

  # update survival steps with rbinom
  pop_next[, 2:nstage] <- rbinom(length(counts), size = counts, prob = surv_vec)

  # update fecundity steps with rpois
  pop_next[, 1] <- rpois(length(fec_vec), lambda = fec_vec)

  # return
  pop_next

}

# internal function: initialise a simulation when inits not provided
initialise <- function(obj, opt, init) {

  if (obj$nspecies > 1) {
    dims <- c(opt$replicates, obj$nclass, obj$nspecies, opt$ntime + 1)
  } else {
    dims <- c(opt$replicates, obj$nclass, opt$ntime + 1)
  }
  ndim <- length(dims)

  # generate initial conditions if not provided
  if (is.null(init)) {

    init <- array(
      options()$aae.pop_initialisation(
        prod(dims[-ndim]), opt$initialise_args
      ),
      dim = dims[-ndim]
    )

  } else {  # check initial values if provided

    # create an error message for re-use
    expected_dims <- dims[1:(ndim - 1)]
    dims_error_msg <- paste0(
      "init has ",
      length(init),
      " elements but must have dimensions (",
      paste0(expected_dims, collapse = ","),
      ") or (",
      paste0(expected_dims[-1], collapse = ","),
      ")"
    )

    # must be a numeric vector or array
    if (!is.numeric(init))
      stop(dims_error_msg, call. = FALSE)

    # if numeric, are the dimensions ok?
    dims_ok <- check_dims(init, expected_dims)

    # error if dims not OK
    if (dims_ok$error)
      stop(dims_error_msg, call. = FALSE)

    # do we need to expand init over replicates?
    if (dims_ok$expand)
      init <- expand_dims(init, expected_dims[1])

  }

  pop <- array(NA, dim = dims)
  pop[seq_len(prod(dims[-ndim]))] <- init

  pop

}

# internal function: check dimensions of initial conditions
check_dims <- function(init, expected_dims) {

  # assume not OK unless inits meet criteria below
  is_ok <- FALSE

  # assume not expanding unless set otherwise
  expand <- FALSE

  # if it's a numeric vector, must have nclass elements
  if (is.null(dim(init))) {

    # are we good?
    if (length(init) == prod(expected_dims[-1])) {
      expand <- TRUE
      is_ok <- TRUE
    }

  } else {

    # are replicates included in init?
    if (length(dim(init)) == length(expected_dims)) {

      if (all.equal(dim(init), expected_dims))
        is_ok <- TRUE

    } else {

      if (length(dim(init)) == (length(expected_dims) - 1)) {
        if (all.equal(dim(init), expected_dims[-1])) {
          expand <- TRUE
          is_ok <- TRUE
        }
      }

    }

  }

  list(error = !is_ok, expand = expand)

}

# internal function: expand initial values
#
#' @importFrom abind abind
expand_dims <- function(init, replicates) {

  # expand over replicates if only one value for each class
  abind::abind(
    lapply(seq_len(replicates), function(x) init),
    along = 0
  )

}

# internal function: initialise simulations with Poisson random draws
#
#' @importFrom stats rpois
initialise_poisson <- function(n, args) {
  do.call(rpois, c(list(n), args))
}

# internal function: set simulation class
as_simulation <- function(x) {
  as_class(x, name = "simulation", type = "array")
}
