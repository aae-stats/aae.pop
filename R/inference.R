#' @name inference
#' @title Perform inference on a population dynamics model
#' @description Estimate parameters in a population dynamics model
#'   defined with the \code{aae.pop} package.
NULL

#' @rdname inference
#'
#' @export
#'
#' @param dynamics a \code{\link{dynamics}} object
#' @param priors a list of priors defined with \code{\link{prior}}
#' @param comparison a function that takes an output from
#'   \code{\link{simulate}} and formats to match inputs
#'   specified in \code{target}
#' @param target values to be compared to simulated model outputs
#'   with \code{comparison}
#' @param method a function call for an inference method that
#'   takes as arguments a function and parameters, returning
#'   a vector of summary statistics that define model fit to data.
#'   Example methods include those defined in the \code{EasyABC}
#'   package
#' @param \dots additional arguments passed to \code{method}
#'
#' @details To be completed.
#'
#' @examples
#' # add
inference <- function(dynamics, priors, comparison, target, method, ...) {

  # collate args into a list
  args <- list(...)

  # define inference_fn
  inference_fn <- define_inference(dynamics, comparison)

  # update parameters
  params <- do.call(method, c(list(inference_fn, priors), args))

  # return
  params

}

#' @rdname inference
#'
#' @export
#'
#' @param dynamics a \code{\link{dynamics}} object
#' @param args a list of arguments for ???
#'
#' @details Return parameters from a dynamics object in format
#'   suitable for inference or for construction of priors.
#'   Flattened parameters include the model matrix in column-major
#'   format, followed by any arguments passed to \code{\link{covariates}},
#'   \code{\link{environmental_stochasticity}},
#'   \code{\link{demographic_stochasticity}},
#'   \code{\link{density_dependence}}, and
#'   \code{\link{density_dependence_n}}, in that order.
#'
#' @examples
#' # add
get_parameters <- function(dynamics, args = list()) {

  # set default arguments passed to dynamic processes
  default_args <- list(
    covariates = list(),
    environmental_stochasticity = list(),
    demographic_stochasticity = list(),
    density_dependence = list(),
    density_dependence_n = list(),
    interaction = list()
  )
  default_args[names(args)] <- args

  # account for multispecies
  if ("multispecies" %in% class(dynamics)) {
    matrix_pars <- lapply(dynamics, function(x) x$dynamics$matrix)
    matrix_pars <- unlist(matrix_pars)
  } else {
    matrix_pars <- dynamics$matrix
    matrix_pars <- c(matrix_pars)
  }

  c(matrix_pars, unlist(default_args))

}

#' @rdname inference
#'
#' @export
#'
#' @param distribution xff
#' @param \dots additional arguments passed to \code{distribution}
#'
#' @details PROBABLY BETTER TO DEFINE UPDATEABLE DEFAULTS
#'   FROM A DYNAMICS OBJECT?
#'
#' @examples
#' # add
prior <- function(distribution, ...) {

  # collate args into a list
  args <- list(...)

  # define and return prior
  do.call(distribution, args)

}

# internal functon: define inference function for dynamics model
define_inference <- function(dynamics, comparison) {

  # return a function that takes in new parameters and
  #   returns summary statistics
  function(params) {

    new_dynamics <- update_dynamics(dynamics, params)

    comparison(new_dynamics)

  }

}

# internal function: update a dynamics object based on a new set of parameters
update_dynamics <- function(dynamics, params) {

  NULL

}
