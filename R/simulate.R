#' @name simulate
#' @title Simulate single or multispecies population dynamics in R
#' @description Simulate population dynamics for one or more
#'   species defined by \code{\link{dynamics}} objects.
NULL

#' @rdname simulate
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
#' @param \dots ignored; included for consistency with \code{simulate} generic
#'   method
#' @param init an array of initial conditions with one row per replicate and one
#'   column per population stage. If \code{obj} has been created with
#'   \code{\link{multispecies}}, initial conditions can be provided as a list or
#'   array with one element or slice per species, or as a matrix, in which case
#'   all species are assumed to share the same initial conditions. Defaults
#'   to \code{NULL}, in which case initial conditions are generated randomly
#'   according to \code{options()$aae.pop_initialisation}
#' @param options a named \code{list} of simulation options. Currently accepted
#'   values are:
#'
#'   - \code{ntime} the number of time steps to simulate, ignored if \code{obj}
#'       includes a \code{\link{covariates}} (default = 50)
#'
#'   - \code{keep_slices} \code{logical} defining whether to keep intermediate
#'       population abundances or (if \code{FALSE}) to return only the final
#'       time slice
#'
#'   - \code{tidy_abundances} a function to handle predicted abundance data
#'       that may be non-integer. Defaults to \code{identity}; suggested
#'       alternatives are \code{floor}, \code{round}, or \code{ceiling}
#'
#'   - \code{initialise_args} a list of arguments passed to the function
#'       used to initialise abundance trajectories. Only used if
#'       \code{init = NULL}. Defaults to \code{options()$aae.pop_lambda},
#'       which specifies lambda for Poisson random draws. The default
#'       initialisation function is defined by
#'       \code{options()$aae.pop_initialisation}.
#'
#'    - \code{update} a function to update abundances from one time
#'       step to the next. Defaults to \code{options()$aae.pop_update}.
#'
#' @param args named list of lists passing arguments to processes defined
#'   in \code{object}, including \code{interaction} for
#'   \code{\link{multispecies}} objects. Lists (up to one per process)
#'   can contain a mix of static, dynamic, and function arguments.
#'   Dynamic arguments must be lists with one element per time step.
#'   Function arguments must be functions that calculate arguments
#'   dynamically in each generation based on from the population dynamics
#'   object, population abundances, and time step in each generation.
#'   All other classes (e.g., single values, matrices, data frames)
#'   are treated as static arguments. Covariates contained in numeric
#'   vectors, matrices, or data frames can be formatted as dynamic
#'   arguments with the \code{format_covariates} function
#' @param args.dyn this argument is deprecated. List of time-varying
#'   values of \code{args}. Defaults to \code{NULL} and requires one
#'   element for each generation (specified by covariates or with
#'   \code{options$ntime}). Dynamic arguments can now be passed directly
#'   to \code{args}
#' @param args.fn this argument is deprecated.
#'   Named list of functions evaluating additional values of
#'   \code{args} based on the matrix and abundances in each generation.
#'   Defaults to \code{NULL}. Function arguments can now be passed
#'   directly to \code{args}
#'
#' @details Includes plot and subset methods
#'
#' @examples
#' # define a population matrix (columns move to rows)
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # can extract standard population matrix summary stats
#' lambda <- Re(eigen(popmat)$values[1])
#'
#' # define a dynamics object
#' dyn <- dynamics(popmat)
#'
#' # simulate from this (50 time steps, 100 replicates)
#' sims <- simulate(dyn, nsim = 100, options = list(ntime = 50))
#'
#' # plot the simulated trajectories
#' plot(sims)
#'
#' # add some density dependence
#' dd <- density_dependence(
#'   masks = reproduction(popmat, dims = 4:5),
#'   funs = ricker(1000)
#' )
#'
#' # update the dynamics object
#' dyn <- update(dyn, dd)
#'
#' # simulate again
#' sims <- simulate(dyn, nsim = 100, options = list(ntime = 50))
#'
#' # and plot
#' plot(sims)
#'
#' # what if we want to add initial conditions?
#' sims <- simulate(
#'   dyn,
#'   init = c(50, 20, 10, 10, 5),
#'   nsim = 100,
#'   options = list(ntime = 50),
#' )
#'
#' # and plot again
#' plot(sims)
#'
#' \dontrun{
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
#'     function(x) rpois(n = length(x), lambda = x),
#'     function(x) rmultiunit(n = 1, mean = x, sd = 0.1 * x)
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
#' # the rmultiunit draws can be slow but we can speed
#' #   this up by calculating them once per generation
#' #   instead of once per replicate within each generation
#' envstoch <- environmental_stochasticity(
#'   masks = list(
#'     reproduction(popmat, dims = 4:5),
#'     transition(popmat)
#'   ),
#'   funs = list(
#'     function(x, ...) rpois(n = length(x), lambda = x),
#'     function(x, mean, sd) {
#'       pnorm(mean + sd * rnorm(length(x)))
#'     }
#'   )
#' )
#'
#' # this requires an argument "function" that takes the
#' #   current state of the population model in each
#' #   iteration and calculates the correct arguments to
#' #   pass to environmental_stochasticty
#' envstoch_function <- function(obj, pop, iter) {
#'   mat <- obj$matrix
#'   if (is.list(mat))
#'     mat <- mat[[iter]]
#'   out <- aae.pop:::unit_to_real(
#'     mat[transition(mat)], 0.1 * mat[transition(mat)]
#'   )
#'   list(mean = out[, 1], sd = out[, 2])
#' }
#'
#' # update the dynamics object and simulate from it
#' dyn <- update(dyn, envstoch)
#' sims <- simulate(
#'   dyn,
#'   init = c(50, 20, 10, 10, 5),
#'   nsim = 100,
#'   args = list(environmental_stochasticity = list(envstoch_function)),
#'   options = list(ntime = 50),
#' )
#'
#' # can also add covariates that influence vital rates
#' #   e.g., a logistic function
#' covars <- covariates(
#'   masks = transition(popmat),
#'   funs = function(mat, x) mat * (1 / (1 + exp(- 10 * x)))
#' )
#'
#' # simulate 50 random covariate values
#' xvals <- matrix(runif(50), ncol = 1)
#'
#' # update the dynamics object and simulate from it.
#' #   Note that ntime is now captured in the 50 values
#' #   of xvals, assuming we pass xvals as an argument
#' #   to the covariates functions
#' dyn <- update(dyn, covars)
#' sims <- simulate(
#'   dyn,
#'   init = c(50, 20, 10, 10, 5),
#'   nsim = 100,
#'   args = list(
#'     covariates = format_covariates(xvals),
#'     environmental_stochasticity = list(envstoch_function)
#'   )
#' )
#'
#' # and can plot these again
#' plot(sims)
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
#'   options = list(update = update_binomial_leslie,
#'                  tidy_abundances = floor),
#'   args = list(
#'     covariates = format_covariates(xvals),
#'     environmental_stochasticity = list(envstoch_function)
#'   )
#' )
#'
#' # and can plot these again
#' plot(sims)
#' }
# nolint start
simulate.dynamics <- function(object,
                              nsim = 1,
                              seed = NULL,
                              ...,
                              init = NULL,
                              options = list(),
                              args = list(),
                              args.dyn = NULL,
                              args.fn = NULL) {
  # nolint end

  # set default options for simulation
  opt <- list(
    ntime = options()$aae.pop_ntime,
    keep_slices = options()$aae.pop_keep_slices,
    tidy_abundances = options()$aae.pop_tidy_abundances,
    initialise_args = list(options()$aae.pop_lambda),
    update = options()$aae.pop_update
  )
  opt[names(options)] <- options

  # check matrix if Leslie updater is used
  if (isTRUE(all.equal(update_binomial_leslie, opt$update))) {
    mat_test <- object$matrix
    mat_test[transition(mat_test)] <- 0
    mat_test[reproduction(mat_test)] <- 0
    if (any(mat_test != 0)) {
      stop("matrix must be a Leslie matrix to use update_binomial_leslie",
           call. = FALSE)
    }
  }

  # set default arguments passed to dynamic processes
  default_args <- list(
    covariates = list(),
    environmental_stochasticity = list(),
    demographic_stochasticity = list(),
    density_dependence = list(),
    density_dependence_n = list(),
    interaction = list()
  )

  # set an identity covariates function if no covariate arguments
  #   provided
  if (is.null(args$covariates))
    object <- use_identity_covariates(object)

  # set identity covariates if args provided and covariates process
  #   is missing
  if (is.multispecies(object)) {

    # do any species have missing covariates?
    no_covariates <- sapply(
      object$dynamics, function(x) is.null(x$covariates)
    )

    # if so, fill with identity covariates
    if (any(no_covariates)) {
      object[no_covariates] <- lapply(
        object[no_covariates],
        use_identity_covariates
      )
    }

  } else {

    if (is.null(object$covariates))
      object$covariates <- use_identity_covariates(object)

  }

  # classify user args by type
  args <- classify_args(args)

  # check args
  if (!is.null(args.dyn) | !is.null(args.fn)) {
    warning("args.dyn and args.fn are deprecated; ",
            "dynamic and function arguments can be ",
            "included directly in args",
            call. = FALSE)
  }

  # set static args here
  default_args[names(args$static)] <- args$static

  # calculate ntime from dynamic args if provided and
  #   check for consistency among dynamic args
  n_dyn_args <- unlist(
    lapply(
      args$dyn,
      function(x) sapply(x, length)
    )
  )
  if (any(n_dyn_args > 0)) {
    n_dyn_args <- n_dyn_args[n_dyn_args > 0]
    if (length(unique(n_dyn_args)) > 1) {
      stop("all dynamic (list) arguments must have the same length",
           call. = FALSE)
    }
    opt$ntime <- unique(n_dyn_args)
  }

  # add nsim into options
  opt$replicates <- nsim

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
  if (is.multispecies(object)) {

    # check if initials are provided as a list with one element
    #   per species
    if (!is.list(init)) {

      # if not, is it an array?
      if (length(dim(init)) == 3) {
        init <- lapply(seq_len(dim(init)[3]), function(i) init[, , i])
      } else { # assume all species shared initial conditions
               #   (possibly NULL)
        init <- lapply(seq_len(object$nspecies), function(i) init)
      }

    }

    # multispecies model, check dims for each
    pop_tmp <- mapply(
      initialise,
      obj = object$dynamics,
      init = init,
      MoreArgs = list(opt = opt, keep_slices = FALSE),
      SIMPLIFY = FALSE
    )

    # do we need to create an object to store everything?
    if (opt$keep_slices) {
      pop <- mapply(
        initialise,
        obj = object$dynamics,
        init = init,
        MoreArgs = list(opt = opt, keep_slices = opt$keep_slices),
        SIMPLIFY = FALSE
      )
    }

  } else {

    # single-species model, check dims
    pop_tmp <- initialise(object, opt, init, keep_slices = FALSE)

    # do we need to create an object to store everything?
    if (opt$keep_slices)
      pop <- initialise(object, opt, init, keep_slices = opt$keep_slices)

  }

  # loop through timesteps, updating population at each timestep
  for (i in seq_len(opt$ntime)) {

    # update args if required
    default_args <- update_args(
      args = args$static,
      dyn = args$dyn,
      fn = args$fn,
      obj = object,
      pop = pop_tmp,
      iter = i
    )

    # split based on multispecies or single species
    if (is.multispecies(object)) {
      pop_tmp <- simulate_once_multispecies(
        iter = i,
        object,
        pop_tmp,
        opt = opt,
        args = default_args
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
        opt = opt,
        args = default_args
      )
      if (opt$keep_slices)
        pop[, , i + 1] <- pop_tmp
    }

  }

  # do we want to keep intermediate abundances or just the final step?
  if (!opt$keep_slices)
    pop <- pop_tmp

  # set appropriate class for outputs
  if (is.multispecies(object)) {
    out <- as_simulation_list(pop)
  } else {
    out <- as_simulation(pop)
  }


  # return
  out

}

#' @importFrom future.apply future_lapply
# internal function: update a single time step for one species
simulate_once <- function(
  iter, obj, pop_t, opt, args, is_expanded = FALSE
) {

  # calculate covariate-altered matrix
  if (is_expanded) {
    mat <- lapply(
      obj$matrix,
      function(x) do.call(
        obj$covariates,
        c(list(x), args$covariates)
      )
    )
  } else {
    mat <- do.call(
      obj$covariates,
      c(list(obj$matrix), args$covariates)
    )
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
        function(x) do.call(
          obj$environmental_stochasticity,
          c(list(x), args$environmental_stochasticity)
        )
      )
    } else {
      mat <- lapply(
        seq_len(opt$replicates),
        function(j) do.call(
          obj$environmental_stochasticity,
          c(list(mat), args$environmental_stochasticity)
        )
      )
      is_expanded <- TRUE
    }
  }

  # tweak matrix to account for density effects on vital rates,
  #   accounting for previously expanded matrix
  if (!is.null(obj[["density_dependence"]])) {
    if (is_expanded) {
      mat <- mapply(
        function(x, y) do.call(
          obj[["density_dependence"]],
          c(list(x, y), args[["density_dependence"]])
        ),
        mat,
        lapply(seq_len(opt$replicates), function(i) pop_t[i, ]),
        SIMPLIFY = FALSE
      )
    } else {
      mat <- lapply(
        seq_len(opt$replicates),
        function(i) do.call(
          obj[["density_dependence"]],
          c(list(mat, pop_t[i, ]), args[["density_dependence"]])
        )
      )
      is_expanded <- TRUE
    }
  }

  # single-step update of abundances
  if (is_expanded) {
    pop_tp1 <- t(mapply(
      opt$update,
      lapply(seq_len(opt$replicates), function(i) pop_t[i, ]),
      mat
    ))
  } else {
    pop_tp1 <- opt$update(pop_t, mat)
  }

  # tweak abundances to add stochastic variation in demographic
  #   outcomes
  if (!is.null(obj$demographic_stochasticity)) {
    pop_tp1 <- t(
      apply(
        pop_tp1,
        1,
        function(x) do.call(
          obj$demographic_stochasticity,
          c(list(x), args$demographic_stochasticity)
        )
      )
    )
  }

  # final hit to abundances if they are rescaled based on biomass
  #   constraints or similar
  if (!is.null(obj$density_dependence_n)) {
    pop_tp1 <- t(
      apply(
        pop_tp1,
        1,
        function(x) do.call(
          obj$density_dependence_n,
          c(list(x), args$density_dependence_n)
        )
      )
    )
  }

  # return tidied abundances (e.g. rounded or floored values)
  opt$tidy_abundances(pop_tp1)

}

# internal function: update a single time step with interacting species
simulate_once_multispecies <- function(iter,
                                       obj,
                                       pop_t,
                                       opt,
                                       args) {

  # vectorised update for all species
  pop_tp1 <- future_lapply(
    seq_len(obj$nspecies),
    simulate_multispecies_internal,
    iter, obj, pop_t, opt, args
  )

  # return tidied abundances (e.g. rounded or floored values)
  opt$tidy_abundances(pop_tp1)

}

# internal function: update one species in a multispecies simulation
#   (to vectorise simulate_once_multispecies)
simulate_multispecies_internal <- function(
  i, iter, obj, pop_t, opt, args
) {

  # pull out relevant object
  dynamics <- obj$dynamics[[i]]

  # and matrix
  mat <- dynamics$matrix

  # rescale matrix according to interspecific interactions
  #   setting a flag to change update step accordingly
  is_expanded <- FALSE
  if (!is.null(obj$interaction[[i]])) {
    mat <- lapply(
      seq_len(opt$replicates),
      function(j)
        do.call(
          obj$interaction[[i]],
          c(list(mat, lapply(pop_t, function(x) x[j, ])), args$interaction)
        )
    )
    is_expanded <- TRUE
  }

  # pass rescaled matrix as part of dynamics
  dynamics$matrix <- mat

  # update and return abundances of species i using single-species updater
  simulate_once(
    iter,
    dynamics,
    pop_t[[i]],
    opt,
    args,
    is_expanded = is_expanded
  )

}

# internal function: update single step of simulation with multiple species
add_multispecies_sims <- function(x, y, iter) {
  x[, , iter + 1] <- y
  x
}

# internal function: initialise a simulation when inits not provided
initialise <- function(obj, opt, init, keep_slices) {

  dims <- c(opt$replicates, obj$nclass, opt$ntime + 1)
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

# internal function: set identity covariates function if covariates
#   are not used
#' @importFrom stats update
use_identity_covariates <- function(obj) {

  # define identity covariates function
  identity_mask <- ifelse(
    is.multispecies(obj),
    obj$dynamics[[1]]$matrix,
    obj$matrix
  )
  identity_covariates <- covariates(
    masks = identity_mask,
    funs = identity
  )

  # loop over species if multispecies, single
  #   update otherwise
  if (is.multispecies(obj)) {
    for (i in seq_len(obj$nspecies)) {
      obj$dynamics <- lapply(
        obj$dynamics,
        update,
        identity_covariates
      )
    }
  } else {
    obj <- update(obj, identity_covariates)
  }

  # return
  obj

}


# internal function: split arguments based on type
classify_args <- function(args) {

  # which arguments are static?
  static <- lapply(args, extract_args, type = "static")

  # which arguments are dynamic?
  dyn <- lapply(args, extract_args, type = "dynamic")

  # which arguments are functions?
  fn <- lapply(args, extract_args, type = "function")

  # and return list of arguments by type
  list(static = static, dyn = dyn, fn = fn)

}

# internal function: extract arguments by type for each
#   process
extract_args <- function(x, type) {

  # work out classes
  arg_class <- sapply(x, class)

  # extract by type:
  #    dynamic if list
  #    function if function,
  #    static otherwise
  if (type == "static")
    x <- x[!arg_class %in% c("list", "function")]
  if (type == "dynamic")
    x <- x[arg_class == "list"]
  if (type == "function")
    x <- x[arg_class == "function"]

  # return
  x

}

# internal function: update arguments based on the current generation
update_args <- function(args, dyn, fn, obj, pop, iter) {

  if (any(sapply(dyn, length) > 0)) {

    # check which exist
    dyn_exist <- names(dyn)

    # update accordingly
    for (i in seq_along(dyn_exist)) {
      for (j in seq_along(dyn[[i]])) {
        args[[dyn_exist[i]]] <- c(
          args[[dyn_exist[i]]], list(dyn[[i]][[j]][[iter]])
        )
      }
    }

  }

  if (any(sapply(fn, length) > 0)) {

    # check which exist
    fn_exist <- names(fn)

    # update accordingly
    for (i in seq_along(fn_exist)) {
      for (j in seq_along(fn[[i]])) {
        fn_eval <- fn[[i]][[j]](obj, pop, iter)
        if (!is.list(fn_eval))
          fn_eval <- list(fn_eval)
        args[[fn_exist[i]]] <- c(args[[fn_exist[i]]], fn_eval)
      }
    }

  }

  # and return
  args

}

# internal function: initialise simulations with Poisson random draws
#
#' @importFrom stats rpois
initialise_poisson <- function(n, args) {
  do.call(rpois, c(list(n), args))
}

# S3 subset method
#' @export
# nolint start
subset.simulation <- function(x, subset, ...) {
  # nolint end
  x <- x[, subset, ]
  as_simulation(x)
}

# S3 subset method
#' @export
# nolint start
subset.simulation_list <- function(x, subset, ...) {
  # nolint end
  for (i in seq_along(x))
    x[[i]] <- x[[i]][, subset, ]
  as_simulation_list(x)
}

# S3 is method
#' @rdname simulate
#'
#' @export
#'
#' @param x an object to pass to \code{is.simulation} or
#'   \code{is.simulation.list}
# nolint start
is.simulation <- function(x) {
  # nolint end
  inherits(x, "simulation")
}

# S3 is method
#' @rdname simulate
#' @export
# nolint start
is.simulation_list <- function(x) {
  # nolint end
  inherits(x, "simulation_list")
}

# S3 print method
#' @export
# nolint start
print.simulation <- function(x, ...) {
  # nolint end
  cat(paste0("Simulated population dynamics for a single species\n"))
}

# S3 print method
#' @export
# nolint start
print.simulation_list <- function(x, ...) {
  # nolint end
  cat(paste0("Simulated population dynamics for ", length(x), " species\n"))
}

# S3 plot method
#' @export
#'
#' @importFrom graphics lines plot
# nolint start
plot.simulation <- function(x, y, ..., class = NULL) {
  # nolint end

  if (is.null(class)) {
    yplot <- apply(x, c(1, 3), sum)
  } else {
    yplot <- x[, class, ]
  }
  xplot <- seq_len(ncol(yplot))
  ylims <- range(yplot)

  args <- list(...)
  plot_defaults <- list(
    x = xplot,
    y = yplot[1, ],
    bty = "l",
    xlab = "Generation",
    ylab = "Abundance",
    las = 1,
    type = "l",
    ylim = ylims
  )
  plot_defaults[names(args)] <- args
  lines_defaults <- list(
    lwd = 1
  )
  lines_defaults[names(args)] <- args

  do.call(plot, plot_defaults)

  for (i in seq_len(nrow(yplot))[-1]) {
    lines_tmp <- c(lines_defaults, list(x = xplot, y = yplot[i, ]))
    do.call(lines, lines_tmp)
  }

}

# S3 plot method
#' @export
# nolint start
plot.simulation_list <- function(x, y, ..., which = seq_along(x)) {
  # nolint end

  nspecies <- length(x)

  for (i in which)
    plot(x[[i]], ...)

}

# internal function: set simulation class
as_simulation <- function(x) {
  type <- "array"
  if (is.matrix(x))
    type <- "matrix"
  as_class(x, name = "simulation", type = type)
}

# internal function: set simulation class
as_simulation_list <- function(x) {

  # each species is a simulation object
  x <- lapply(x, as_simulation)

  # but combination of species is a list
  as_class(x, name = "simulation_list", type = "list")

}
