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

  as_mask(row(matrix) == 1 & col(matrix) != 1 & col(matrix) %in% dims)

}

#' @rdname masks
#'
#' @export
survival <- function(matrix, dims = NULL) {

  if (is.null(dims))
    dims <- seq_len(ncol(matrix))

  as_mask(row(matrix) == col(matrix) & col(matrix) %in% dims)

}

#' @rdname masks
#'
#' @export
transition <- function(matrix, dims = NULL) {

  if (is.null(dims))
    dims <- seq_len(ncol(matrix) - 1)

  as_mask(
    row(matrix) == (col(matrix) + 1) & col(matrix) %in% dims
  )

}

#' @rdname masks
#'
#' @export
all_cells <- function(matrix, dims = NULL) {
  as_mask(matrix(TRUE, nrow = nrow(matrix), ncol = ncol(matrix)))
}

#' @rdname masks
#'
#' @export
all_stages <- function(matrix, dims = NULL) {

  nstage <- ncol(matrix)

  if (is.null(dims))
    dims <- seq_len(nstage)

  as_mask(seq_len(nstage) %in% dims)

}

#' @rdname masks
#'
#' @export
#'
#' @param \dots sj
combine <- function(...) {
  UseMethod("combine")
}

#' @rdname masks
#'
#' @export
#'
#' @details something
#'
#' @examples
#' # add
combine.mask <- function(...) {
  dots <- list(...)
  masks <- abind::abind(dots, along = 3)
  apply(masks, c(1, 2), any)
}

#' @rdname masks
#'
#' @export
combine.function <- function(...) {
  dots <- list(...)
  function(matrix, dims = NULL) {
    out <- list()
    for (i in seq_along(dots))
      out[[i]] <- dots[[i]](matrix, dims)
    combine(out)
  }
}

#' @rdname masks
#'
#' @export
combine.default <- function(...) {
  dots <- list(...)
  classes <- sapply(dots, class)
  stop("combine is not defined for objects of class",
       clean_paste(classes),
       call. = FALSE)
}

# internal function: set mask class
as_mask <- function(x) {
  as_class(x, name = "mask", type = "matrix")
}
