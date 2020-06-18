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
#' @param args named list of lists passing arguments to processes defined
#'   in \code{object}.
#' @param \dots currently ignored, included for consistency with simulate
#'   generic
#'
#' @details To be completed.
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
                              options = list(),
                              args = list(),
                              ...) {
  # nolint end

  # set default options for simulation
  opt <- list(
    ntime = options()$aae.pop_ntime,
    keep_slices = options()$aae.pop_keep_slices,
    tidy_abundances = options()$aae.pop_tidy_abundances,
    initialise_args = list(options()$aae.pop_lambda)
  )
  opt[names(options)] <- options

  # set default arguments passed to dynamic processes
  default_args <- list(
    covariates = list(),
    environmental_stochasticity = list(),
    demographic_stochasticity = list(),
    density_dependence = list(),
    density_dependence_n = list()
  )

  # add nsim into options
  opt$replicates <- nsim

  # use the number of covariate values instead of fixed ntime if
  #   covariates are provided
  if (object$nspecies > 1) {
    if (object$include_covariates)
      opt$ntime <- object$ntime
  } else {
    if (!is.null(object$covariates))
      opt$ntime <- object$ntime
  }

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
  if (object$nspecies > 1) {
    pop_tmp <- lapply(object$dynamics, initialise, opt, init, keep_slices = FALSE)
    if (opt$keep_slices)
      pop <- lapply(object$dynamics, initialise, opt, init, keep_slices = opt$keep_slices)
  } else {
    pop_tmp <- initialise(object, opt, init, keep_slices = FALSE)
    if (opt$keep_slices)
      pop <- initialise(object, opt, init, keep_slices = opt$keep_slices)
  }

  # loop through timesteps, updating population at each timestep
  for (i in seq_len(opt$ntime)) {

    # split based on multispecies or single species
    if (object$nspecies > 1) {
      pop_tmp <- simulate_once_multispecies(
        iter = i,
        object,
        pop_tmp,
        opt = opt
      )
      if (opt$keep_slices) {
        pop <- mapply(
          add_multispecies_sims, pop, pop_tmp,
          MoreArgs = list(iter = i),
          SIMPLIFY = FALSE
        )
      }
    } else {
      pop_tmp <- simulate_once(
        iter = i,
        object,
        pop_tmp,
        opt = opt
      )
      if (opt$keep_slices)
        pop[, , i + 1] <- pop_tmp
    }

  }

  # do we want to keep intermediate abundances or just the final step?
  if (!opt$keep_slices)
    pop <- pop_tmp

  # set appropriate class for outputs
  if (object$nspecies > 1) {
    out <- as_simulation_list(pop)
  } else {
    out <- as_simulation(pop)
  }


  # return
  out

}

#' @importFrom future.apply future_lapply future_mapply
# internal function: update a single time step for one species
simulate_once <- function(iter, obj, pop_t, opt, is_expanded = FALSE, args = list()) {

  # matrix will be a list if expanded over covariates
  if (!is.null(obj$covariates)) {
    mat <- obj$matrix[[iter]]
  } else {
    mat <- obj$matrix
  }

  # keep pop_t as a matrix if replicates == 1
  if (opt$replicates == 1)
    pop_t <- matrix(pop_t, nrow = 1)

  # draw stochastic matrix values if env stoch included,
  #   accounting for previously expanded matrix if multispecies
  if (!is.null(obj$environmental_stochasticity)) {
    if (is_expanded) {
      mat <- lapply(
        mat,
        function(x) do.call(obj$environmental_stochasticity, c(list(x), args$environmental_stochasticity))
      )
    } else {
      mat <- lapply(
        seq_len(opt$replicates),
        function(j) do.call(obj$environmental_stochasticity, c(list(mat), args$environmental_stochasticity))
      )
      is_expanded <- TRUE
    }
  }

  # tweak matrix to account for density effects on vital rates,
  #   accounting for previously expanded matrix
  if (!is.null(obj$density_dependence)) {
    if (is_expanded) {
      mat <- mapply(
        obj$density_dependence,
        mat,
        lapply(seq_len(opt$replicates), function(i) pop_t[i, ]),
        SIMPLIFY = FALSE
      )
    } else {
      mat <- lapply(
        seq_len(opt$replicates),
        function(i) obj$density_dependence(mat, pop_t[i, ])
      )
      is_expanded <- TRUE
    }
  }

  # single-step update of abundances
  if (is_expanded) {
    pop_tp1 <- t(mapply(
      options()$aae.pop_update,
      lapply(seq_len(opt$replicates), function(i) pop_t[i, ]),
      mat
    ))
  } else {
    pop_tp1 <- options()$aae.pop_update(pop_t, mat)
  }

  # tweak abundances to add stochastic variation in demographic
  #   outcomes
  if (!is.null(obj$demographic_stochasticity))
    pop_tp1 <- t(apply(pop_tp1, 1, function(x) do.call(obj$demographic_stochasticity, c(list(x), args$demographic_stochasticity))))

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
                                       opt) {

  # vectorised update for all species
  pop_tp1 <- future_lapply(
    seq_len(obj$nspecies),
    simulate_multispecies_internal,
    iter, obj, pop_t, opt
  )

  # return tidied abundances (e.g. rounded or floored values)
  opt$tidy_abundances(pop_tp1)

}

# internal function: update one species in a multispecies simulation
#   (to vectorise simulate_once_multispecies)
simulate_multispecies_internal <- function(i, iter, obj, pop_t, opt) {

  # pull out relevant object
  dynamics <- obj$dynamics[[i]]

  # matrix will be a list if expanded over covariates
  if (!is.null(dynamics$covariates)) {
    mat <- dynamics$matrix[[iter]]
  } else {
    mat <- dynamics$matrix
  }

  # rescale matrix according to interspecific interactions
  #   setting a flag to change update step accordingly
  is_expanded <- FALSE
  if (!is.null(obj$interaction[[i]])) {
    if (is_expanded) {
      mat <- mapply(
        obj$interaction[[i]],
        mat,
        lapply(seq_len(opt$replicates),
               function(j) lapply(pop_t, function(x) x[j, ])),
        SIMPLIFY = FALSE
      )
    } else {
      mat <- lapply(
        seq_len(opt$replicates),
        function(j) obj$interaction[[i]](mat, lapply(pop_t, function(x) x[j, ]))
      )
      is_expanded <- TRUE
    }
  }

  # include correct matrix in dynamics and set covariates to NULL
  #   because it's already handled above
  dynamics$matrix <- mat
  dynamics$covariates <- NULL

  # update and return abundances of species i using single-species updater
  simulate_once(iter, dynamics, pop_t[[i]], opt, is_expanded = is_expanded)

}

# internal function: update single step of simulation with multiple species
add_multispecies_sims <- function(x, y, iter) {
  x[, , iter + 1] <- y
  x
}

# internal function: initialise a simulation when inits not provided
initialise <- function(obj, opt, init, keep_slices) {

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

  # only return single slice (initials) if !keep_slices
  if (!keep_slices)
    pop <- array(pop[seq_len(prod(dims[-ndim]))], dim = dims[-ndim])

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

# internal function: set simulation class
as_simulation_list <- function(x) {

  # each species is a simulation object
  x <- lapply(x, as_simulation)

  # but combination of species is a list
  as_class(x, name = "simulation_list", type = "list")

}
