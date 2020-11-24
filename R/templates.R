#' @name templates
#' @title Parameterised population dynamics objects
#' @description Use pre-defined population dynamics objects to
#'   define a matrix model for a species with known parameters.
NULL

# internal method
get_template <- function(sp, ...) {

  # unpack dots
  arg_list <- list(...)

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)
  all_parameters <- do.call(get(paste0("template_", sp)), arg_list)

  # return collated dynamics object
  do.call(dynamics, all_parameters)

}

#' @rdname templates
#'
#' @export
#'
#' @importFrom stats rnorm
#'
#' @param k carrying capacity
#' @param reproductive integer vector specifying reproductive
#'   age classes or stages in the population matrix
#' @param system ecosystem type defining elements of population
#'   dynamics. Currently implemented for \code{macquarie_perch},
#'   which has different covariate effects in lakes and rivers
#'
#' @param \dots additional arguments passed to templates or args
#'   functions
#'
#' @details These functions (e.g. \code{murray_cod}) return a collated
#'   \code{\link{dynamics}} object parameterised with values based
#'   on existing data sets and published works. Currently implemented
#'   species are: Murray cod (*Maccullochella peelii*) and Macquarie
#'   perch (*Macquaria australasica*).
#'
#'   In some cases, additional arguments might need to be passed to
#'   \code{\link{simulate}}, and it may be useful to define these as
#'   part of the template. The \code{get_args} function
#'   handles this situation. Arguments should be
#'   specified with \code{args_} followed by
#'   a species name or identifier (e.g \code{args_my_species}).
#'   Arguments functions should return a series of named lists
#'   for any of \code{args}, \code{args.dyn}, or \code{args.fn} (see
#'   \code{\link{simulate}} for descriptions of these terms).
#'
#' @examples
#' # define a basic model for Murray cod with
#' #   carrying capacity = 25000
#' mc <- murray_cod(k = 25000)
#'
#' # simulate from this model
#' sims <- simulate(mc, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
murray_cod <- function(k = 20000, ...) {
  get_template(sp = "murraycod", k = k, ...)
}

#' @rdname templates
#'
#' @export
#'
#' @importFrom stats pnorm rnorm runif
macquarie_perch <- function(
  k = 1000,
  reproductive = 3:30,
  system = "lake",
  ...
) {
  get_template(
    sp = "macquarieperch",
    k = k,
    reproductive = reproductive,
    system = system,
    ...)
}

# internal function: define species defaults
template_murraycod <- function(k = 20000) {

  # how many stages are we going to work with?
  nstage <- 25

  # define base matrix
  mat <- matrix(0, nrow = nstage, ncol = nstage)
  survival_mask <- combine(
    transition(mat), survival(mat, dims = nstage)
  )
  mat[survival_mask] <- c(
    0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
    0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
    0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253, 0.9375
  )
  yoy_surv <- 0.5 * 0.0122 * 0.1225
  reproduction_mask <- reproduction(mat, dims = c(5:nstage))
  mat[reproduction_mask] <- yoy_surv * c(
      3000, 5000, 7000, 9000, 12000, 16000,
      20000, 25000, 30000, 34000, 38000,
      41000, 43000, 45000, 47000, 48000,
      48000, 49000, 49000, 49000, 50000
    )

  # define basic BH density dependence
  biomass_dd <- function(k, dims) {
    function(x, n) {
      sum_n <- sum(n[min(dims):length(n)])
      x * ifelse(sum_n > k, k / sum_n, 1)
    }
  }
  dd_stages <- list(
    c(3:4),
    c(5:7),
    c(8:10),
    c(11:14),
    c(15:25)
  )
  biomass_dd_list <- lapply(dd_stages, biomass_dd, k = k)
  biomass_mask_list <- lapply(
    dd_stages, transition, mat = mat
  )
  biomass_mask_list[[length(biomass_mask_list)]][nstage, nstage] <- TRUE
  dd_fns <- c(
    list(beverton_holt(k = k)), biomass_dd_list
  )
  dd_masks <- c(
    list(reproduction(mat, dims = c(5:nstage))),
    biomass_mask_list
  )
  dd <- density_dependence(dd_masks, dd_fns)

  # basic single variable covariate function
  cov_funs <- function(mat, x) {
    mat * (1 / (1 + exp(-0.5 * (x + 10))))
  }
  cov_masks <- transition(mat)
  covars <- covariates(
    masks = cov_masks,
    funs = cov_funs
  )

  # define environmental stochasticity based on known standard deviations of
  #   parameters
  # survival
  survival_env <- function(x) {

    # define base parameters for SD
    sd_scale <- c(0.2, rep(0.15, 3), rep(0.1, 3), rep(0.075, 3), rep(0.05, 15))
    sd_tmp <- sd_scale * c(
      0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
      0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
      0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253, 0.9375
    )

    # simulate normal random variates
    out <- rnorm(length(x), mean = x, sd = sd_tmp)

    # check none sit outside 0/1
    out <- ifelse(out > 1, 1, out)
    out <- ifelse(out < 0, 0, out)

    # return
    out

  }
  # reproduction
  reproduction_env <- function(x) {

    # define base parameters for SD
    yoy_surv <- 0.5 * 0.0122 * 0.1225
    sd_tmp <- yoy_surv * c(
      1500, 2400, 3200, 4000, 5300, 6900,
      8400, 10500, 12600, 13900, 15600,
      16800, 17600, 18500, 18800, 19200,
      19200, 19600, 19600, 19600, 20000
    )

    # simulate normal random variates
    out <- rnorm(length(x), mean = x, sd = sd_tmp)

    # check none are negative (but can be > 1)
    out <- ifelse(out < 0, 0, out)

    # return
    out

  }
  # collate into a enviro stochasticity object
  envstoch <- environmental_stochasticity(
    masks = list(survival_mask, reproduction_mask),
    funs = list(survival_env, reproduction_env)
  )

  # set demographic stochasticity
  # survival
  demo_fn <- function(x) {
    rpois(length(x), lambda = x)
  }
  demostoch <- demographic_stochasticity(
    masks = all_classes(mat),
    funs = demo_fn
  )

  # return template
  list(
    matrix = mat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = demostoch,
    density_dependence = dd,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_troutcod <- function() {

  # return template
  list(
    matrix = NULL,
    covariates = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_goldenperch <- function() {

  # return template
  list(
    matrix = NULL,
    covariates = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_silverperch <- function() {

  # return template
  list(
    matrix = NULL,
    covariates = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_macquarieperch <- function(
  k = 1000,                     # adult female carrying capacity
  reproductive = 3:30,          # reproductive age classes
  system = "lake"               # define covariate type by system
) {

  # set default system
  if (length(system) > 1)
    system <- system[1]

  # and check it
  if (!system %in% c("lake", "river"))
    stop("system must be one of lake or river", call. = FALSE)

  # define a survival function, adding dots to soak up extra
  #    arguments within simulate
  survival_gen <- function(
    mat, mean_real, sd_real, perfect_correlation = TRUE, ...
  ) {
    rmultiunit_from_real(
      n = 1,
      mean_real = mean_real,
      sd_real = sd_real,
      perfect_correlation = perfect_correlation
    )
  }

  # define a reproduction function, being careful with argument
  #   names to avoid conflicts with any arguments in
  #   survival_gen, which would then require multiple different
  #   arguments with the same name in simulate
  reproduction_gen <- function(
    mat,
    fec_mean = c(1.68, -0.302, 2.886),
    fec_sd = c(0.3, 0.05, 0.15),
    early_mean,
    early_sd,
    recruit_failure = 0.25,
    contributing_min = 0.5,
    contributing_max = 1.0,
    ...
  ) {

    # generate stochastic values for early life
    #   survival (eggs, larvae, young-of-year)
    early_surv <- rmultiunit_from_real(n = 1, mean = early_mean, sd = early_sd)

    # otherwise draw random variates for the three model parameters
    y1 <- rnorm(n = 1, mean = fec_mean[1], sd = fec_sd[1])
    y2 <- rnorm(n = 1, mean = fec_mean[2], sd = fec_sd[2])
    y3 <- rnorm(n = 1, mean = fec_mean[3], sd = fec_sd[3])

    # generate reproduction estimates for all adult age classes, incorporating
    #   stochastic early life estimates
    y2_term <- exp(y2 %o% reproductive)
    y1_y2 <- log(
      43.15 * exp(sweep(y2_term, 1, -y1, "*"))
    )
    reprod <- exp(sweep(2.295 * y1_y2, 1, y3, "+"))

    # did recruitment fail?
    recruit_binary <- ifelse(recruit_failure >= runif(1), 0, 1)

    # make contributing a stochastic variables
    if (contributing_min > contributing_max) {
      contributing_min <- contributing_max
      warning(
        "contributing_min was greater than contributing_max; ",
        "both values have been set to contributing_max",
        call. = FALSE
      )
    }
    contributing <- runif(1, min = contributing_min, max = contributing_max)

    # add early life survival and muliply by 0.5
    #   to account for a 50:50 sex ratio
    0.5 * recruit_binary * contributing * reprod * prod(early_surv)

  }

  # define mean survival
  survival_mean <- c(
    0.25, 0.44, 0.56, 0.63, 0.69, 0.72, 0.75, 0.78, 0.79,
    0.81, 0.82, 0.83, 0.83, 0.84, 0.84, 0.84, 0.85, 0.85, 0.84,
    0.84, 0.84, 0.83, 0.82, 0.80, 0.78, 0.76, 0.71, 0.63, 0.48
  )

  # calculate reproduction at mean parameter values
  reproduction_mean <- exp(
    2.295 *
      log(43.15 * exp(- 1.68 * exp(-0.302 * reproductive))) +
      2.886) *
    0.5 *                        # 50:50 sex ratio
    prod(c(0.5, 0.013, 0.13))    # add early life survival

  # define population matrix
  nclass <- length(survival_mean) + 1
  popmat <- matrix(0, nrow = nclass, ncol = nclass)
  popmat[transition(popmat)] <- survival_mean
  popmat[reproduction(popmat, dims = reproductive)] <- reproduction_mean

  # define density dependence, only affects adult survival
  #   and reproductive stages
  density_masks <- list(
    transition(popmat, dims = reproductive),
    reproduction(popmat, dims = reproductive)
  )

  # top-down effects of competition for habitat, at
  #   carrying capacity k
  topdown_fn <- function(mat, pop, ...) {
    sum_n <- sum(pop[reproductive])
    ifelse(sum_n > k, k / sum_n, 1) * mat
  }

  # positive density dependence (Allee effect)
  allee_fn <- function(mat, pop, allee_strength = 1, allee_factor = 10, ...) {
    sum_n <- sum(pop[reproductive])
    allee <- (2 / (1 + exp(-sum_n / (allee_strength * allee_factor)))) - 1
    allee * mat
  }

  # and collate masks and functions in a single object
  dens_depend <- density_dependence(
    masks = density_masks,
    funs = list(topdown_fn, allee_fn)
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      transition(popmat),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(survival_gen, reproduction_gen)
  )

  # use density_dependence_n to include stocking or
  #   translocations (removals)
  dd_n_masks <- list(
    all_classes(popmat, dim = 1),
    all_classes(popmat, dim = 2),
    all_classes(popmat, dim = reproductive)
  )
  dd_n_fns <- list(
    function(pop, n_yoy, ...)
      add_remove(pop = pop, n = n_yoy, add = TRUE),
    function(pop, n_twoplus, ...)
      add_remove(pop = pop, n = n_twoplus, add = TRUE),
    function(pop, n_adult, ...)
      add_remove(pop = pop, n = n_adult, add = TRUE)
  )
  dens_depend_n <- density_dependence_n(
    masks = dd_n_masks,
    funs = dd_n_fns
  )

  # define covariate effects on recruitment in all systems
  recruit_effects <- function(mat, x, ...) {

    # effect of spawning-flow magnitude on recruitment
    # negative effect of Nov/Dec discharge on recruitment
    log_flow <- log(x$spawning_flow + 0.01)
    scale_factor <- exp(-0.1 * log_flow - 0.1 * (log_flow ^ 2))
    scale_factor[scale_factor > 1] <- 1
    scale_factor[scale_factor < 0] <- 0
    mat <- mat * scale_factor

    # effect of spawning-flow variability on recruitment
    # negative effect of Nov/Dec discharge variability on recruitment
    #   (days with more than 100% change from previous)
    mat <- mat * exp(-0.05 * x$spawning_variability)

    # return
    mat

  }

  # define covariate effects on adult survival in all systems
  survival_effects <- NULL

  # define covariate masks in all systems
  recruit_masks <- reproduction(popmat)
  survival_masks <- NULL

  # add system-specific covariate effects
  if (system == "lake") {

    # negative effect of rising lake level on YOY
    recruit_effects_lake <- function(mat, x, ...) {
      mat * (1 / (1 + exp(-0.5 * (x$water_level_change + 10))))
    }

    # update effects and masks
    recruit_effects <- list(recruit_effects, recruit_effects_lake)
    recruit_masks <- list(recruit_masks, reproduction(popmat))

  }
  if (system == "river") {

    # negative effect of rising river level on YOY
    recruit_effects_river <- function(mat, x, ...) {
      mat * (1 / (1 + exp(-0.5 * (x$river_height_change + 10))))
    }

    # negative effect of low flows on adult survival
    survival_effects_river <- function(mat, x, ...) {

      # positive effect of flow on overall population growth rate
      #    (based on individuals >= 1 year old)
      log_flow <- log(x$average_daily_flow + 0.01)
      scale_factor <- exp(0.3 * log_flow - 0.3 * (log_flow ^ 2))
      scale_factor[scale_factor > 1] <- 1
      scale_factor[scale_factor < 0] <- 0
      mat <- mat * scale_factor

      # return
      mat

    }

    # update effects and masks
    recruit_effects <- list(recruit_effects, recruit_effects_river)
    recruit_masks <- list(recruit_masks, reproduction(popmat))
    survival_effects <- list(survival_effects, survival_effects_river)
    survival_masks <- list(
      survival_masks, transition(popmat, dims = reproductive)
    )

  }

  # compile covariates process
  covars <- covariates(
    masks = list(recruit_masks, survival_masks),
    funs = list(recruit_effects, survival_effects)
  )

  # nolint start
  # print a message to ensure args are included
  message('Arguments are required to simulate Macquarie perch dynamics.\n',
          'Default arguments can be accessed with get_args("macquarie_perch").')
  # nolint end

  # return
  list(
    matrix = popmat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = NULL,
    density_dependence = dens_depend,
    density_dependence_n = dens_depend_n
  )

}

# internal function: define species defaults
template_australiansmelt <- function() {

  # return template
  list(
    matrix = NULL,
    covariates = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_commongalaxias <- function() {

  # return template
  list(
    matrix = NULL,
    covariates = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

#' @rdname templates
#'
#' @export
get_args <- function(sp, ...) {

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)

  # initialise NULL output
  out <- NULL

  # check if species has arguments available
  available <- exists(
    paste0("args_", sp),
    envir = getNamespace("aae.pop"),
    mode = "function"
  )

  # collapse dots into a list
  arg_list <- list(...)

  # get arguments if available
  if (available)
    out <- do.call(get(paste0("args_", sp)), arg_list)

  # return
  out

}

# internal function: define macquarie perch arguments
args_macquarieperch <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  allee_strength = 1,
  contributing_min = 0.75,
  contributing_max = 1.0,
  recruit_failure = 0
) {

  # expand n, start, end if required
  if (length(n) == 1)
    n <- rep(n, 3)
  if (length(start) == 1)
    start <- rep(start, 3)
  if (length(end) == 1)
    end <- rep(end, 3)

  # check for other lengths of n, start, or end
  if (any(c(length(n), length(start), length(end)) != 3)) {
    stop("n, start, and end must be vectors with 1 or 3 elements",
         call. = FALSE)
  }

  # helper to define removals process
  define_removals <- function(start, end, n, add = add) {

    # set up a sequence of iterations at which individuals are removed
    if (!is.list(n)) {
      n <- mapply(
        zeros_and_fill,
        n,
        start,
        end,
        MoreArgs = list(len = ntime),
        SIMPLIFY = FALSE
      )
    } else {
      if (!all.equal(sapply(n, length), ntime)) {
        stop("if n is a list, each element must be a vector ",
             "with one value for each time step",
             call. = FALSE)
      }
    }

    # define this as a function
    translocate <- function(obj, pop, iter) {

      # return
      list(
        n_yoy = n[[1]][iter],
        n_twoplus = n[[2]][iter],
        n_adult = n[[3]][iter],
        add = add
      )

    }

    # return
    translocate

  }

  # helper to calculate real-valued parameters for survival
  #   simulation
  transform_survival <- function(obj, pop, iter) {

    # pull out the population matrix in the current time step
    mat <- obj$matrix
    if (is.list(mat))
      mat <- mat[[iter]]

    # wrap up all survival means and SDs, including early life
    #  (this allows a single call to `unit_to_real`, which is slow)
    survival_mean <- c(
      0.5, 0.013, 0.13,     # early life
      mat[transition(mat)]  # from population matrix in current time step
    )
    survival_sd <- c(
      0.1, 0.007, 0.028,  # early life
      0.05, 0.09, 0.11, 0.10, 0.10, 0.07, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.06, 0.05
    )

    # convert unit interval to real line equivalents
    out <- unit_to_real(
      unit_mean = survival_mean,
      unit_sd = survival_sd
    )

    # separate early life from other estimates
    idx <- seq_len(nrow(out)) > 3

    # return
    list(
      mean_real = out[idx, 1],    # for survival_gen
      sd_real = out[idx, 2],      # for survival_gen
      early_mean = out[!idx, 1],  # for reproduction_gen
      early_sd = out[!idx, 2]     # for reproduction_gen
    )

  }

  # define static elements of each scenario
  args_list <- list(

    # set as 1 (default) or 2
    density_dependence = list(allee_strength = allee_strength),

    # set contributing as random uniform on 0.75-1.0 by default
    # set recruit_failure at 0 by default
    environmental_stochasticity = list(
      contributing_min = contributing_min,
      contributing_max = contributing_max,
      recruit_failure = recruit_failure
    )

  )

  # define functions to update simulate args dynamically
  args_function <- list(

    # to pre-transform unit to real and back
    environmental_stochasticity = transform_survival,

    # to include additions or removals of individuals
    density_dependence_n = define_removals(
      start = start, end = end, n = n, add = add
    )
  )

  # return named list of args
  list(
    static = args_list,
    funs = args_function
  )

}

# internal function to initialise a vector of zeros
#   and fill a subset of values
zeros_and_fill <- function(x, start, end, len) {
  out <- rep(0, len)
  out[start:end] <- x
  out
}

# internal function to handle translocations or stocking
add_remove <- function(
  pop,
  n,
  add = TRUE,
  ...
) {

  # only add/remove individuals if required
  if (n > 0) {

    # check there are enough individuals to remove
    if (sum(pop) < n & !add) {
      n <- sum(pop)
      warning("removal required more individuals than ",
              "were available; reduced to ",
              sum(pop),
              " individuals",
              call. = FALSE)
    }

    # are we removing?
    if (!add) {

      # if so, expand n to remove from random age classes
      n_by_age <- rep(seq_along(pop), times = pop)
      idx <- sample.int(
        length(n_by_age), size = n, replace = FALSE
      )
      n_change <- table(n_by_age[idx])

    } else {

      # otherwise, add to age classes at random
      n_change <- table(
        sample(seq_along(pop), size = n, replace = TRUE)
      )

    }

    n_expanded <- rep(0, length(pop))
    names(n_expanded) <- as.character(seq_along(pop))
    n_expanded[names(n_change)] <- n_change

    # are we adding or removing?
    if (add)
      n_expanded <- -n_expanded

    # update pop abundances
    pop <- pop - n_expanded

  }

  # return
  pop

}

# convert species names to standardised common names
#  (could replace with sci names if ambiguous)
parse_species <- function(sp) {

  # create lookup table with common synonyms
  sp_list <- list(
    "murraycod" = "murraycod",
    "maccullochellapeelii" = "murraycod",
    "troutcod" = "troutcod",
    "maccullochellamacquariensis" = "troutcod",
    "goldenperch" = "goldenperch",
    "macquariaambigua" = "goldenperch",
    "silverperch" = "silverperch",
    "bidyanusbidyanus" = "silverperch",
    "macquarieperch" = "macquarieperch",
    "macquariaaustralasica" = "macquarieperch",
    "australiansmelt" = "australiansmelt",
    "retropinnasemoni" = "australiansmelt",
    "commongalaxias" = "commongalaxias",
    "galaxiasmaculatus" = "commongalaxias"
  )

  # clean sp to remove underscores or spaces or hyphens
  sp_clean <- gsub("_|-| ", "", sp)

  # pull out a match if there, otherwise check the template
  #   function exists and return as is (error otherwise)
  if (sp_clean %in% names(sp_list)) {
    sp <- sp_list[sp_clean]
  } else {
    if (!exists(paste0("template_", sp))) {
      stop(sp, " does not have a matching template function",
           call. = FALSE)
    }
  }

  # return appropriate species
  sp

}
