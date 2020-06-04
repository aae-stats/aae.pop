#' @name stochasticity
#' @title Specify environmental and demographic stochasticity in models of
#'   population dynamics
#' @description Specify environmental stochasticity (random variation
#'   in vital rates) and demographic stochasticity (random variation
#'   in population outcomes).
NULL

#' @rdname stochasticity
#'
#' @export
#'
#' @param masks dfds
#' @param funs d
#'
#' @details something
#'
#' @examples
#' # add
environmental_stochasticity <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(x) {
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]])
      x
    }
  } else {
    fn <- function(x) {
      do_mask(x, masks, funs)
    }
  }

  as_environmental_stochasticity(fn)

}

#' @rdname stochasticity
#'
#' @export
#'
#' @details something
#'
#' @examples
#' # add
demographic_stochasticity <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(x) {
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]])
      x
    }
  } else {
    fn <- function(x) {
      do_mask(x, masks, funs)
    }
  }

  as_demographic_stochasticity(fn)

}

# internal function: set environmental_stochasticity class
as_environmental_stochasticity <- function(x) {
  as_class(
    x, name = "environmental_stochasticity", type = "function"
  )
}

# internal function: set demographic_stochasticity class
as_demographic_stochasticity <- function(x) {
  as_class(
    x, name = "demographic_stochasticity", type = "function"
  )
}
