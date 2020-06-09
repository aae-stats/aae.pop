#' @name multispecies
#' @title Create a population dynamics object with multiple species
#' @description Define population dynamics for multiple species from
#'   a set of single-species \code{\link{dynamics}} objects and
#'   defined pairwise interactions.
NULL

#' @rdname multispecies
#'
#' @export
#'
#' @param \dots \code{pairwise_interaction} objects defining all
#'   pairwise interactions of interest
multispecies <- function(...) {

  # collate dots into list
  dots <- list(...)

  # all ... should be pairwise interactions objects
  classes <- sapply(dots, function(x) class(x)[1])
  if (!all(classes == "interaction")) {
    stop("all inputs to multispecies should be interaction ",
         "objects created with pairwise_interaction",
         call. = FALSE)
  }

  # pull out dynamics$hex for all and use to find unique objects
  dynamics <- get_unique_dynamics(dots)

  # return
  as_multispecies(
    list(dynamics = dynamics)
  )

}

#' @rdname multispecies
#'
#' @export
#'
#' @param target population whose vital rates are affected by the
#'   pairwise interaction
#' @param source population whose abundances affect the vital rates
#'   of \code{target}
#' @param masks masks defining which vital rates are influenced by
#'   each function
#' @param funs functions that take vital rates and abundances of the
#'   \code{source} population as inputs and return scaled vital rates
#'
#' @details DKLDJF
pairwise_interaction <- function(target, source, masks, funs) {

  # define function that specifies effects of source on target
  if (is.list(masks)) {
    effects <- function(x, n) {
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]], n)
      x
    }
  } else {
    effects <- function(x, n) {
      do_mask(x, masks, funs, n)
    }
  }

  # return
  as_interaction(
    list(target = target,
         source = source,
         effects = effects)
  )

}

# internal function: extract unique dynamics objects from
#   list of interaction objects
get_unique_dynamics <- function(interactions) {

  # which hex values do we have?
  hex_target <- sapply(interactions, function(x) x$target$hex)
  hex_source <- spaply(interactions, function(x) x$source$hex)

  # create a vector of unique targets and sources
  hex_list <- unique(c(hex_target, hex_source))

  # which pops match these unique hexes?
  target_ids <- match(hex_list, hex_target)
  out <- lapply(interactons[target_ids], function(x) x$target)
  if (any(is.na(target_ids))) {
    source_ids <- match(hex_list[is.na(target_ids)], source_ids)
    out <- c(
      out,
      lapply(interactions[source_ids], function(x) x$source)
    )
  }

  # return all required dynamics objects
  out

}

# internal function: set multispecies class
as_multispecies <- function(x) {

  # inherit from dynamics
  x <- as_dynamics(x)

  # and then add multispecies class on top of that
  as_class(x, name = "multispecies", type = "dynamics")

}

# internal function: set interaction class
as_interaction <- function(x) {
  as_class(x, name = "interaction", type = "list")
}
