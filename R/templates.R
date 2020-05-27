#' @name templates
#' @title Use a fully parameterised population dynamics object
#' @description Use pre-defined population dynamics objects to
#'   define a matrix model for a species with known parameters,
#'   optionally including new covariates.
NULL

#' @rdname templates
#'
#' @export
#'
#' @param sp string specifying species common name (see details for
#'   a list of species currently included)
#' @param x a vector or matrix of covariates used to define
#'   time-varying vital rates according to a pre-defined
#'   functional response
#' @param \dots additional objects passed to \code{\link{dynamics}},
#'   overwriting the default, pre-defined values for these classes.
#'   Must be one or more of \code{\link{covariates}},
#'   \code{\link{environmental_stochasticity}},
#'   \code{\link{demographic_stochasticity}},
#'   \code{\link{density_dependence}}, or
#'   \code{\link{density_dependence_n}}
#'
#' @details The \code{get_template} function and associated
#'   wrapper functions (e.g. \code{murraycod}) return a collated
#'   \code{\link{dynamics}} object parameterised with values based
#'   on existing data sets and published works. Currently implemented
#'   species are: Murray cod (*Maccullochella peelii*),
#'   trout cod (*Maccullochella macquariensis*), golden perch
#'   (*Macquaria ambigua*), silver perch (*Bidyanus bidyanus*),
#'   Macquarie perch (*Macquaria australasica*), Australian smelt
#'   (*Retropinna semonia*), and common galaxias (*Galaxias maculatus*).
#'
#'   The \code{get_template} function can be used with user-defined
#'   templates. Templates must be functions with a name defined by
#'   \code{template_} followed by a species name or identifier (e.g.
#'   \code{template_my_species}). These templates must return a named
#'   list with at least a population matrix model (\code{matrix}), and
#'   other arguments to pass to \code{\link{dynamics}}. Terms
#'   specifying \code{covariates} would typically be defined
#'   in a template with the function only (\code{covariate_function}),
#'   so that templates can be used with different covariate sets.
#'   If defined in this way, \code{get_template("my_species")} will
#'   return a compiled dynamics object for my_species.
#'
#' @examples
#' # add
get_template <- function(sp, x = NULL, ...) {

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)
  all_parameters <- get(paste0("template_", sp))()

  # optional: add covariates
  if (!is.null(x)) {
    all_parameters$covariates <-
      covariates(x = x, fun = all_parameters$covariate_function)
  }

  # remove covariate_function from all_parameters in all cases
  all_parameters$covariate_function <- NULL

  # unpack dots and replace defaults if any objects provided
  arg_list <- list(...)
  if (length(arg_list) > 0) {
    arg_types <- sapply(arg_list, function(x) class(x)[1])
    all_parameters[arg_types] <- arg_list
  }

  # return collated dynamics object
  do.call(dynamics, all_parameters)

}

#' @rdname templates
#'
#' @export
murraycod <- function(x = NULL, ...) {
  get_template(sp = "murraycod", x = x, ...)
}

# internal function: define species defaults
template_murraycod <- function() {

  # how many stages are we going to work with?
  nstage <- 25

  # define base matrix
  mat <- matrix(0, nrow = nstage, ncol = nstage)
  survival_mask <- combine_masks(
    transition(mat), col(mat) == row(mat) & col(mat) == nstage
  )
  mat[survival_mask] <- c(
    0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
    0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
    0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253, 0.9375
  )
  mat[nstage, nstage] <- 0.9375
  yoy_surv <- 0.5 * 0.0122 * 0.1225
  reproduction_mask <- reproduction(mat, dims = c(5:nstage))
  mat[reproduction_mask] <- yoy_surv * c(
      3000, 5000, 7000, 9000, 12000, 16000,
      20000, 25000, 30000, 34000, 38000,
      41000, 43000, 45000, 47000, 48000,
      48000, 49000, 49000, 49000, 50000
    )

  # define basic BH density dependence
  dd <- beverton_holt(
    masks = reproduction(mat, dims = c(5:nstage)),
    params = list(K = 20000)
  )

  # work out more complex, biomass-based density dependence
  #  (clustered if/else statements)
  # mask is survival, function is
  # ifelse(K < summed_abunds, K / summed_abunds, 1)
  # Variation is in class clusters, e.g. 3+4, 5-7, 8-10, 11-14, 15-25. Clusters
  #   affect all things younger but nothing older
  #   (e.g. 3+4 = sum(abunds[3:nstage], 8-10 = sum(abunds[8:nstage]))

  # basic single variable covariate function
  covs <- function(mat, x) {
    mat[transition(mat)] <- mat[transition(mat)] * (1 / (1 + exp(-0.5 * x)))
  }

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
    masks = combine_masks(survival_mask, reproduction_mask),
    funs = demo_fn
  )

  # return template
  list(
    matrix = mat,
    covariate_function = covs,
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
    covariate_function = NULL,
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
    covariate_function = NULL,
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
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_macquarieperch <- function() {

  # return template
  list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

}

# internal function: define species defaults
template_australiansmelt <- function() {

  # return template
  list(
    matrix = NULL,
    covariate_function = NULL,
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
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
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
