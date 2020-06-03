#' @name dispersal
#' @title Specify dispersal between populations in a metapopulation model
#' @description Specify disperal between populations, including stochasticity
#'   and density dependence in dispersal parameters
NULL

#' @rdname dispersal
#'
#' @export
#'
#' @param stochasticity_masks dfds
#' @param stochasticity_funs dfds
#' @param density_masks dfds
#' @param density_funs dfds
#'
#' @details something
#'
#' @examples
#' # add
dispersal <- function(stochasticity_masks,
                      stochasticity_funs,
                      density_masks,
                      density_funs) {

  if (is.list(stochasticity_masks)) {
    stoch_fn <- function(x) {
      for (i in seq_along(stochasticity_masks)) {
        x[stochasticity_masks[[i]]] <-
          stochasticity_funs[[i]](x[stochasticity_masks[[i]]])
      }
      x
    }
  } else {
    stoch_fn <- function(x) {
      x[stochasticity_masks] <-
        stochasticity_funs(x[stochasticity_masks])
      x
    }
  }

  if (is.list(density_masks)) {
    dens_fn <- function(x) {
      for (i in seq_along(density_masks)) {
        x[density_masks[[i]]] <-
          density_funs[[i]](x[density_masks[[i]]])
      }
      x
    }
  } else {
    dens_fn <- function(x) {
      x[density_masks] <-
        density_funs(x[density_masks])
      x
    }
  }

  fn_list <- list(
    stochasticity = stoch_fn,
    density = dens_fn
  )

  as_dispersal(fn_list)

}

# internal function: set demographic_stochasticity class
as_dispersal <- function(x) {
  as_class(x, name = "dispersal", type = "list")
}

