#' @name density_functions
#' @title Common forms of density dependence
#' @description Use pre-defined forms of density dependence
#'   based on common density-dependence functions.
NULL

#' @rdname density_functions
#'
#' @export
#'
#' @param k carrying capacity used to define models of
#'   density dependence. See details for
#'   currently implemented models and their parameters.
#' @param exclude vector of classes to exclude from calculation
#'   of total population density. Defaults to NULL, in which
#'   case all classes are used
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence. Currently implemented models
#'   are the Ricker model and Beverton-Holt model, both with
#'   a single parameter \code{k}.
#'
#' @examples
#' # add
beverton_holt <- function(k, exclude = NULL) {

  function(x, n) {
    if (!is.null(exclude))
      n <- n[-exclude]
    x / (1 + x * sum(n) / k)
  }

}

#' @rdname density_functions
#'
#' @export
#'
#' @examples
#' # add
ricker <- function(k, exclude = NULL) {

  function(x, n) {
    if (!is.null(exclude))
      n <- n[-exclude]
    x * exp(1 - sum(n) / k) / exp(1)
  }

}
