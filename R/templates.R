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
macquarie_perch <- function(k = 1000, ...) {
  get_template(sp = "macquarieperch", k = k, ...)
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
template_macquarieperch <- function(k = 1000) {

  # need some helper functions within template
  #
  # simulate a vector of perfectly correlated
  #   survival values with known mean and sd
  survival_gen <- function(n, mean, sd, perfect_correlation = TRUE) {

    # how many parameters are we dealing with?
    npar <- length(mean)

    # simulate random values from a standard normal
    z_variates <- matrix(rnorm(npar * n), ncol = n)

    # do we want perfectly correlated or uncorrelated?
    if (perfect_correlation) {
      out <- t(pnorm(mean + sd %o% z_variates[1, ]))
    } else {
      out <- t(pnorm(mean + sweep(z_variates, 1, sd, "*")))
    }

    # return
    out

  }

  # simulate a vector of perfectly correlated
  #   fecundity values with known mean and sd
  fecundity <- function(age,
                        mean,
                        egg_survival,
                        larval_survival,
                        yoy_survival,
                        n = 1,
                        sd = NULL,
                        recruit_failure = 0.25,
                        contributing_min = 0.5,
                        contributing_max = 1.0) {

    # is SD provided?
    if (is.null(sd)) {

      # if not, we want to return the mean fecundity only
      y1 <- mean$y1
      y2 <- mean$y2
      y3 <- mean$y3

    } else {

      # otherwise draw random variates for the three model parameters
      y1 <- rnorm(n = n, mean = mean$y1, sd = sd$y1)
      y2 <- rnorm(n = n, mean = mean$y2, sd = sd$y2)
      y3 <- rnorm(n = n, mean = mean$y3, sd = sd$y3)

    }

    # calculate fecundity
    y2_term <- exp(y2 %o% age)
    y1_y2 <- log(
      43.15 * exp(sweep(y2_term, 1, -y1, "*"))
    )
    fec <- exp(sweep(2.295 * y1_y2, 1, y3, "+"))

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

    # now add in the remaining terms and return
    fec *
      contributing *
      recruit_binary *
      egg_survival *
      larval_survival *
      yoy_survival *
      0.5

  }

  # set main parameters (nclass)
  nclass <- 30

  # define survival parameters
  survival_params <- list(
    mean = c(0.13, 0.25, 0.44, 0.56, 0.63, 0.69, 0.72, 0.75, 0.78, 0.79,
             0.81, 0.82, 0.83, 0.83, 0.84, 0.84, 0.84, 0.85, 0.85, 0.84,
             0.84, 0.84, 0.83, 0.82, 0.80, 0.78, 0.76, 0.71, 0.63, 0.48),
    sd = c(0.028, 0.05, 0.09, 0.11, 0.10, 0.10, 0.07, 0.08, 0.08, 0.08,
           0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
           0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.06, 0.05),
    egg_mean = 0.5,
    egg_sd = 0.1,
    larvae_mean = 0.013,
    larvae_sd = 0.007
  )

  # define parameters for fecundity
  fecundity_params <- list(
    mean = list(y1 = 1.68, y2 = -0.302, y3 = 2.886),
    sd = list(y1 = 0.3, y2 = 0.05, y3 = 0.15)
  )

  # define fecundity
  fec_mean <- fecundity(
    age = 4:30,
    mean = fecundity_params$mean,
    egg_survival = survival_params$egg_mean,
    larval_survival = survival_params$larvae_mean,
    yoy_survival = survival_params$mean[1],
    recruit_failure = 0,
    contributing_min = 1.0,
    contributing_max = 1.0
  )

  # define population matrix
  popmat <- matrix(0, nrow = nclass, ncol = nclass)
  popmat[transition(popmat, dims = 1:29)] <- survival_params$mean[2:30]
  popmat[reproduction(popmat, dims = 4:30)] <- fec_mean

  # define density dependence
  dens_masks <- list(
    transition(popmat, dims = 4:29),
    reproduction(popmat, dims = 4:30)
  )
  density_fn <- function(mat, pop, ...) {
    sum_n <- sum(pop[4:30])
    ifelse(sum_n > k, k / sum_n, 1) * mat
  }
  allee_fn <- function(mat, pop, allee_strength = 1, allee_factor = 10, ...) {
    sum_n <- sum(pop[4:30])
    allee <- (2 / (1 + exp(-sum_n / (allee_strength * allee_factor)))) - 1
    mat <- allee * mat
    mat
  }
  dens_fns <- list(
    density_fn,
    allee_fn
  )
  dens_depend <- density_dependence(
    masks = dens_masks,
    funs = dens_fns
  )

  # define environmental stochasticity
  fec_envstoch <- function(mat,
                           egg_mean,
                           egg_sd,
                           larvae_mean,
                           larvae_sd,
                           yoy_mean,
                           yoy_sd,
                           contributing_min = 0.5,
                           contributing_max = 1.0,
                           recruit_failure = 0.25,
                           ...) {

    young_surv <- survival_gen(n = 1,
                               mean = c(egg_mean, larvae_mean, yoy_mean),
                               sd = c(egg_sd, larvae_sd, yoy_sd))
    fecundity(
      age = 4:30,
      egg_survival = young_surv[1],
      larval_survival = young_surv[2],
      yoy_survival = young_surv[3],
      mean = fecundity_params$mean,
      sd = fecundity_params$sd,
      recruit_failure = recruit_failure,
      contributing_min = contributing_min,
      contributing_max = contributing_max
    )
  }
  envstoch_masks <- list(
    transition(popmat, dims = 1:29),
    reproduction(popmat, dims = 4:30)
  )
  envstoch_fns <- list(
    function(mat, mean, sd, ...) survival_gen(n = 1, mean = mean, sd = sd),
    fec_envstoch
  )
  envstoch <- environmental_stochasticity(
    masks = envstoch_masks,
    funs = envstoch_fns
  )

  # add density_dependence_n to deal with translocations
  macperch_dens_depend_n <- function(
    pop,
    translocate = FALSE,
    n_translocate = NULL,
    add = TRUE,
    adults = 4:30
  ) {

    # only remove individuals if required
    if (translocate) {

      # set to zero removed if not specified
      if (is.null(n_translocate))
        n_translocate <- rep(0, 4)

      # check there are enough individuals to translocate (if removing)
      if (sum(pop[adults]) < n_translocate[4] & !add) {
        n_translocate[4] <- sum(pop[adults])
        warning("translocation required more adults than ",
                "were available. Reducing the number of ",
                "translocated individuals to ",
                sum(pop[adults]),
                call. = FALSE)
      }

      # are we removing?
      if (!add) {

        # if so, expand n to remove from random age classes
        n_adult_by_age <- rep(adults, times = pop[adults])
        adult_idx <- sample.int(
          length(n_adult_by_age), size = n_translocate[4], replace = FALSE
        )
        n_adult <- table(n_adult_by_age[adult_idx])

      } else {

        # otherwise, add to age classes at random
        n_adult <- table(
          sample(adults, size = n_translocate[4], replace = TRUE)
        )

      }
      n_expanded <- rep(0, length(adults))
      names(n_expanded) <- as.character(adults)
      n_expanded[names(n_adult)] <- n_adult
      n_take <- c(n_translocate[1],
                  n_translocate[2],
                  n_translocate[3],
                  n_expanded)

      # are we adding or removing?
      if (add)
        n_take <- -n_take

      # check that there are enough YOY and 1+/2+ individuals to take
      if (n_take[1] > pop[1]) {
        n_take[1] <- pop[1]
        warning("translocation required more YOY than ",
                "were available; ",
                "translocation reduced to ",
                pop[1],
                " individuals",
                call. = FALSE)
      }
      if (n_take[2] > pop[2]) {
        n_take[2] <- pop[2]
        warning("translocation required more 1+ individuals than ",
                "were available; ",
                "translocation reduced to ",
                pop[2],
                " individuals",
                call. = FALSE)
      }
      if (n_take[3] > pop[3]) {
        n_take[3] <- pop[3]
        warning("translocation required more 2+ individuals than ",
                "were available; ",
                "translocation reduced to ",
                pop[3],
                " individuals",
                call. = FALSE)
      }

      # update pop abundances
      pop <- pop - n_take

    }

    # return
    pop

  }
  dens_depend_n <- density_dependence_n(masks = all_classes(popmat),
                                        funs = macperch_dens_depend_n)


  # define covariate effects
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

    # negative effect of lake level and increasing
    #   temperature on growth of YOY, plus a
    #   negative effect of increasing (rising) dam height
    mat <- mat * (1 / (1 + exp(-0.5 * (x$water_level_change + 10))))

    # return
    mat

  }
  adult_effects <- function(mat, x, ...) {

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
  covar_funs <- list(
    recruit_effects,
    adult_effects
  )
  covar_masks <- list(
    reproduction(popmat),
    transition(popmat, dims = 4:30)
  )
  covars <- covariates(
    masks = covar_masks,
    funs = covar_funs
  )

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
  n = c(0, 0, 0, 0),
  start = 1,
  end = 1,
  add = TRUE,
  allee_strength = 1,
  contributing_min = 0.75,
  contributing_max = 1.0,
  recruit_failure = 0) {

  # helper to define translocation process
  define_translocation <- function(start, end, n_translocate, add = TRUE) {

    # set up a sequence of iterations at which individuals are removed
    iter_seq <- start:end

    # define this as a function
    translocate <- function(obj, pop, iter) {

      # initialise
      translocate <- FALSE

      # flag to pass to dens_depend_n
      if (iter %in% iter_seq & sum(n_translocate) > 0)
        translocate <- TRUE

      # return
      list(
        translocate = translocate,
        n_translocate = n_translocate,
        add = add
      )

    }

    # return
    translocate

  }

  # define an args function so we can pass transformed values to envstoch
  transform_survival <- function(obj, pop, iter) {

    sd_vals <- c(
      0.1, 0.007,
      0.028, 0.05, 0.09, 0.11, 0.10, 0.10, 0.07, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.06, 0.05
    )
    mat <- obj$matrix
    if (is.list(mat))
      mat <- mat[[iter]]

    mean_vals <- c(0.5,     # mean egg survival
                   0.013,   # mean larval survival
                   0.13,    # mean YOY survival
                   mat[transition(mat, dims = 1:29)])

    out <- unit_to_real(
      unit_mean = mean_vals,
      unit_sd = sd_vals
    )

    nmat <- nrow(out)

    # return
    list(mean = out[4:nmat, 1],
         sd = out[4:nmat, 2],
         egg_mean = out[1, 1],
         egg_sd = out[1, 2],
         larvae_mean = out[2, 1],
         larvae_sd = out[2, 2],
         yoy_mean = out[3, 1],
         yoy_sd = out[3, 2]
    )

  }

  # define static elements of each scenario
  args_list <- list(

    # set as 1 (default) or 2
    density_dependence = list(allee_strength = allee_strength),

    # set contributing as random uniform on 0.5-1.0 by default
    # set recruit_failure as 0.25 (default), 0.3, or 0.5
    environmental_stochasticity = list(
      contributing_min = contributing_min,
      contributing_max = contributing_max,
      recruit_failure = recruit_failure
    )

  )

  # define functions to update simulate args dynamically
  #   - environmental stochasticity: necessary to speed up
  #     transformation of survival from unit to real and back
  #   - density_dependence_n: necessary to include translocation
  args_function <- list(

    # to pre-transform unit to real and back
    environmental_stochasticity = transform_survival,

    # options: 100-600 YOY for 5 or 10 years
    #          50-125 1+ for 5 or 10 years
    #          50-100 2+ per year
    #          5 or 15 adults for 1 or 5 years
    density_dependence_n = define_translocation(
      start = start, end = end, n = n, add = add
    )
  )

  # return named list of args
  list(
    static = args_list,
    funs = args_function
  )


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
