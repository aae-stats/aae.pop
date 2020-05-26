#' @name density_dependence
#' @title Specify density dependence in models of population dynamics
#' @description Specify density dependence in vital rates
#'   (\code{density_dependence}) and in total abundances
#'   (\code{density_dependence_n}).
NULL

#' @rdname density_dependence
#'
#' @export
#'
#' @param masks s
#' @param funs dfd
#'
#' @details something
#'
#' @examples
#' # add
density_dependence <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(x, n) {
      for (i in seq_along(masks))
        x[masks[[i]]] <- funs[[i]](x[masks[[i]]], n)
      x
    }
  } else {
    fn <- function(x, n) {
        x[masks] <- funs(x[masks], n)
        x
    }
  }

  as_density_dependence(fn)

}

#' @rdname density_dependence
#'
#' @export
#'
#' @details something
#'
#' @examples
#' # add
density_dependence_n <- function(masks, funs) {

  if (is.list(masks)) {
    fn <- function(pop_t) {
      for (i in seq_along(masks))
        pop_t[masks[[i]]] <- funs[[i]](pop_t[masks[[i]]])
      pop_t
    }
  } else {
    fn <- function(pop_t) {
      pop_t[masks] <- funs(pop_t[masks])
      pop_t
    }
  }

  as_density_dependence_n(fn)

}

#' @rdname density_dependence
#'
#' @export
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence, such as the Ricker model and
#'   the Beverton-Holt model.
#'
#' @examples
#' # add
ricker <- function(masks) {
  density_dependence(masks, lapply(masks, ricker_dd))
}

#' @rdname density_dependence
#'
#' @export
#'
#' @param params K
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence, such as the Ricker model and
#'   the Beverton-Holt model.
#'
#' @examples
#' # add
beverton_holt <- function(masks, params) {

  if (is.list(masks)) {
    out <- density_dependence(masks, lapply(masks, function() bh_dd(params)))
  } else {
    out <- density_dependence(masks, bh_dd(params))
  }

  out

}

# internal function: define ricker density dependence
ricker_dd <- function(x) {
  NULL
}

# internal function: define Beverton-Holt density dependence
bh_dd <- function(params) {
  function(x, n) {
    x / (1 + x * sum(n) / params$K)
  }
}

# internal function: set density_dependence class
as_density_dependence <- function(x) {
  as_class(x, name = "density_dependence", type = "function")
}

# internal function: set density_dependence_n class
as_density_dependence_n <- function(x) {
  as_class(x, name = "density_dependence_n", type = "function")
}
