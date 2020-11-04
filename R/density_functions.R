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
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence. Currently implemented models
#'   are the Ricker model and Beverton-Holt model, both with
#'   a single parameter \code{k}.
#'
#' @examples
#' # add
beverton_holt <- function(k) {

  function(x, n) {
    x / (1 + x * sum(n) / k)
  }

}

#' @rdname density_functions
#'
#' @export
#'
#' @examples
#' # add
ricker <- function(k) {

  function(x, n) {
    x * exp(1 - sum(n) / k) / exp(1)
  }

}

