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
#' @param matrix a population dynamics matrix for which
#'   a particular mask is required. Only used to determine
#'   mask dimensions, so can be any matrix with appropriate
#'   dimensions
#' @param dims a numeric value or vector identifying
#'   subsets of cells to include in a given mask
#'
#' @details To be completed.
#'
#' @examples
#' # define a population
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # pull out reproductive elements
#' reproduction(popmat)
#'
#' # what if only 4 and 5 year olds reproduce?
#' reproduction(popmat, dims = 4:5)
#'
#' # define survival elements
#' survival(popmat)
#'
#' # what if 1 and 2 year olds always transition?
#' survival(popmat, dims = 3:5)
#'
#' # and transitions
#' transition(popmat)
#'
#' # combine transitions and reproduction of 4 and 5 year olds
#' combine(reproduction(popmat, dims = 4:5), transition(popmat))
#'
#' # can also mask the population vector in this way
#' # pull out all classes
#' all_classes(popmat)
#'
#' # and just 3-5 year olds
#' all_classes(popmat, dims = 3:5)
reproduction <- function(matrix, dims = NULL) {
  if (is.null(dims)) {
    dims <- seq_len(ncol(matrix))
  }

  as_mask(row(matrix) == 1 & col(matrix) != 1 & col(matrix) %in% dims)
}

#' @rdname masks
#'
#' @export
survival <- function(matrix, dims = NULL) {
  if (is.null(dims)) {
    dims <- seq_len(ncol(matrix))
  }

  as_mask(row(matrix) == col(matrix) & col(matrix) %in% dims)
}

#' @rdname masks
#'
#' @export
transition <- function(matrix, dims = NULL) {
  if (is.null(dims)) {
    dims <- seq_len(ncol(matrix) - 1)
  }

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
all_classes <- function(matrix, dims = NULL) {
  nclass <- ncol(matrix)

  if (is.null(dims)) {
    dims <- seq_len(nclass)
  }

  as_mask(
    matrix(seq_len(nclass) %in% dims, ncol = 1)
  )
}

#' @rdname masks
#'
#' @export
#'
#' @param \dots a set of masks or masking functions
#'   to be combined into a single mask by one of the
#'   \code{combine} methods
combine <- function(...) {
  UseMethod("combine")
}

# S3 method
#' @export
combine.mask <- function(...) {
  # turn dots into a list
  dots <- list(...)

  # check classes of all dots
  classes <- sapply(dots, function(x) class(x)[1])

  # error if classes not OK
  if (!all(classes == "mask")) {
    stop("combine is not defined for mask objects combined with ",
      "objects of class ",
      clean_paste(classes[!(classes %in% c("mask"))]),
      call. = FALSE
    )
  }

  # return combined mask if all OK
  masks <- abind::abind(dots, along = 3)
  as_mask(apply(masks, c(1, 2), any))
}

# S3 method
#' @export
combine.function <- function(...) {
  # turn dots into a list
  dots <- list(...)

  # check classes of all dots
  classes <- sapply(dots, function(x) class(x)[1])

  # error if classes not OK
  if (!all(classes == "function")) {
    stop("combine is not defined for function objects combined with ",
      "objects of class ",
      clean_paste(classes[!(classes %in% c("function"))]),
      call. = FALSE
    )
  }

  # return function if all OK
  function(matrix, dims = NULL) {
    out <- list()
    for (i in seq_along(dots)) {
      out[[i]] <- dots[[i]](matrix, dims)
    }
    do.call(combine, out)
  }
}

# S3 method
#' @export
combine.default <- function(...) {
  # turn dots into a list
  dots <- list(...)

  # check classes of all dots
  classes <- sapply(dots, function(x) class(x)[1])

  # error and tell user which class was passed
  stop("combine is not defined for objects of class ",
    clean_paste(classes[!(classes %in% c("mask", "function"))]),
    call. = FALSE
  )
}

# internal function: set mask class
as_mask <- function(x) {
  as_class(x, name = "mask", type = "matrix")
}
