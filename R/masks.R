#' @name masks
#' @title Isolate elements of population dynamics models
#' @description Helper functions to isolate particular components
#'   of a population dynamics model, such as the reproduction terms,
#'   transition/growth terms, or particular life stages from an
#'   abundance vector, such as pre- or post-reproductive stages.
NULL

#' @rdname masks
#'
#' @export
#'
#' @param matrix s
#' @param dims dfd
#'
#' @details something
#'
#' @examples
#' # add
reproduction <- function(matrix, dims = NULL) {

  if (is.null(dims))
    dims <- seq_len(ncol(matrix))

  row(matrix) == 1 & col(matrix) != 1 & col(matrix) %in% dims

}

#' @rdname masks
#'
#' @export
survival <- function(matrix, dims = NULL) {

  if (is.null(dims))
    dims <- seq_len(ncol(matrix))

  row(matrix) == col(matrix) & col(matrix) %in% dims

}

#' @rdname masks
#'
#' @export
transition <- function(matrix, dims = NULL) {

  if (is.null(dims))
    dims <- seq_len(ncol(matrix) - 1)

  row(matrix) == (col(matrix) + 1) & col(matrix) %in% dims

}

#' @rdname masks
#'
#' @export
all_cells <- function(matrix) {
  matrix(TRUE, nrow = nrow(matrix), ncol = ncol(matrix))
}

#' @rdname masks
#'
#' @export
all_stages <- function(matrix, dims = NULL) {

  nstage <- ncol(matrix)

  if (is.null(dims))
    dims <- seq_len(nstage)

  seq_len(nstage) %in% dims

}

#' @rdname masks
#'
#' @export
#'
#' @details something
#'
#' @examples
#' # add
combine_masks <- function(...) {
  dots <- list(...)
  masks <- abind::abind(dots, along = 3)
  apply(masks, c(1, 2), any)
}

# internal function: set mask class
as_mask <- function(x) {
  as_class(x, name = "mask", type = "matrix")
}
