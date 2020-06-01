#' @name metapopulation
#' @title Create a population dynamics object with multiple populations
#' @description Define population dynamics for multiple populations of
#'   a single species linked by dispersal.
NULL

#' @rdname metapopulation
#'
#' @export
metapopulation <- function(structure, dynamics, dispersal, ...) {

  # structure should be a matrix with nothing on diagonal and 0/1 T/F elsewhere
  npop <- nrow(structure)

  # dynamics can be a list of length npop or a single matrix,
  #   expand if single matrix

  # work out underlying stage info
  nstage <- nrow(dynamics[[1]]$matrix)

  # could define indexing that works with npop and nstage to give each
  #   element of structure in expanded dims (needs rows + cols)

  # create banded diagonal with dynamics on the band

  # add dispersal to off-diagonal elements where structure == 1

  # set outputs from input dynamics objects
  #  - decide how to handle conflicting settings for input pops
  out <- NULL
  as_metapopulation(out)

}

# internal function: set metapopulation class
as_metapopulation <- function(x) {
  as_class(x, name = "metapopulation", type = "list")
}
