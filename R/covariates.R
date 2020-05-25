#' @name covariates
#'
#' @title Specify covariate dependence in models of population dynamics
#'
#' @description Specify relationship between a vector or matrix of covariates
#'   and vital rates.
#'
#' @export
#'
#' @param x fk
#' @param fun fk
#'
#' @details something
#'
#' @examples
#' # add
covariates <- function(x, fun) {

  # don't precalculate this step -- means same modifier can be used for
  #    several matrices with different dimensions (assuming fun is general enough for this, e.g.,
  #    defined with masks not explicit dims)

  # convert all x to matrix with time slices in columns
  if (!is.matrix(x))
    x <- matrix(x, ncol = 1)

  # return
  as_covariates(list(x = x, fun = fun, ntime = nrow(x)))

}

# internal function: set covariates class
as_covariates <- function(x) {
  as_class(x, name = "covariates", type = "list")
}

