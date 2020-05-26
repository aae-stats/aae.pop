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
#' @examples
#' # add
get_template <- function(sp, x = NULL, ...) {

  # draw up relevant parameters
  all_parameters <- templates(sp)

  # optional: add covariates
  if (!is.null(x))
    all_parameters$covariates <- covariates(x = x, fun = all_parameters$covariate_function)

  # remove covariate_function from all_parameters in all cases
  all_parameters$covariate_function <- NULL

  # unpack dots and replace defaults with any provided objects
  arg_list <- list(...)
  arg_types <- sapply(arg_list, function(x) class(x)[1])
  all_parameters[arg_types] <- arg_list

  # return collated dynamics object
  do.call(dynamics, all_parameters)

}

#' @rdname templates
#'
#' @export
murraycod <- function(x = NULL, ...) {
  get_template(sp = "murraycod", x = x, ...)
}

# internal information for dynamics templates
templates <- function(sp) {

  # define templates for each species

  # Murray cod
  nstage <- 25
  mc_matrix <- matrix(0, nrow = nstage, ncol = nstage)
  mc_matrix[transition(mc_matrix)] <- c(
    0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
    0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
    0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253
  )
  mc_matrix[nstage, nstage] <- 0.9375
  survival_to_yoy <- 0.5 * 0.0122 * 0.1225
  mc_matrix[reproduction(mc_matrix, dims = c(5:nstage))] <-
    survival_to_yoy * c(
    3000, 5000, 7000, 9000, 12000, 16000, 20000, 25000, 30000, 34000,
    38000, 41000, 43000, 45000, 47000, 48000, 48000, 49000, 49000,
    49000, 50000
  )
  mc_dd <- beverton_holt(
    masks = reproduction(mc_matrix, dims = c(5:nstage)),
    params = list(K = 20000)
  )
  ## PICK UP HERE AND WORK IT OUT
  mc_covs <- function(mat, x) {
    mat[transition(mat)] <- mat[transition(mat)] * (1 / (1 + exp(-5 * x)))
  }
  murraycod <- list(
    matrix = murraycod_matrix,
    covariate_function = mc_covs,
    environmental_stochasticity = NULL, # define from SDs of parameters
    demographic_stochasticity = NULL,   # think about ways to use rbin()
    density_dependence = mc_dd,
    density_dependence_n = NULL
  )

  # trout cod
  troutcod <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # golden perch
  goldenperch <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # silver perch
  silverperch <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # Macquarie perch
  macquarieperch <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # Australian smelt
  australiansmelt <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # common galaxias
  commongalaxias <- list(
    matrix = NULL,
    covariate_function = NULL,
    environmental_stochasticity = NULL,
    demographic_stochasticity = NULL,
    density_dependence = NULL,
    density_dependence_n = NULL
  )

  # collate all into a single list
  template_list <- list(
    "murraycod" = murraycod,
    "troutcod" = troutcod,
    "goldenperch" = goldenperch,
    "silverperch" = silverperch,
    "macquarieperch" = macquarieperch,
    "australiansmelt" = australiansmelt,
    "commongalaxias" = commongalaxias
  )

  # return
  template_list[sp]

}


