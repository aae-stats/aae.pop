# package file

#' aae.pop: flexible simulation of multispecies population dynamics
#' @name aae.pop
#'
#' @description aae.pop supports flexible specification of population
#'   dynamics with tools for simulating single and multispecies dynamics,
#'   including in metapopulation structures.
#'
#'   See \code{\link{dynamics}} and \code{\link{simulate}} for detailed
#'   examples of model definition and simulation.
#'
#' @keywords internal
"_PACKAGE"

# set some options when package is loaded
.onLoad <- function(libname, pkgname) { # nolint

  # set update type for single-step updates
  options(aae.pop_update = update_crossprod)

  # set default initialisation method
  options(
    aae.pop_initialisation = initialise_poisson,
    aae.pop_lambda = 10
  )

  # set default options for fecundity classes
  reproduction_default <- function(x) {
    reproduction(x, dims = 2:ncol(x))
  }
  options(aae.pop_reproduction_mask = reproduction_default)

  # set default options for simulations
  options(
    aae.pop_ntime = 50,
    aae.pop_keep_slices = TRUE,
    aae.pop_tidy_abundances = identity
  )
}
