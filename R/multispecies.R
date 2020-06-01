#' @name multispecies
#' @title Create a population dynamics object with multiple species
#' @description Define population dynamics for multiple species from
#'   a set of single-species \code{\link{dynamics}} objects and
#'   defined pairwise interactions.
NULL

#' @rdname multispecies
#'
#' @export
multispecies <- function(interactions, ...) {

  # all ... should be dynamics objects


  # pull out dynamics$hex for all and use to find unique objects

  # set outputs from input dynamics objects
  #  - decide how to handle conflicting settings for input pops
  out <- NULL

  as_multispecies(out)

}

#' @rdname multispecies
#'
#' @export
pairwise_interaction <- function(target, cause, masks, funs) {

  out <- NULL

  as_interaction(out)

}

# internal function: set multispecies class
as_multispecies <- function(x) {
  as_class(x, name = "multispecies", type = "list")
}

# internal function: set interaction class
as_interaction <- function(x) {
  as_class(x, name = "interaction", type = "list")
}
