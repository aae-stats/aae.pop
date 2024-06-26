#' @name multispecies
#' @title Create a population dynamics object with multiple species
#' @description Define population dynamics for multiple species from
#'   a set of single-species \code{\link{dynamics}} objects and
#'   defined pairwise interactions.
#'
#' @export
#'
#' @param \dots \code{pairwise_interaction} objects defining
#'   a set of pairwise interactions between species
multispecies <- function(...) {

  # collate dots into list
  dots <- list(...)

  # all ... should be pairwise interactions objects
  classes <- sapply(dots, function(x) class(x)[1])
  if (!all(classes == "interaction")) {
    stop("all inputs to multispecies should be interaction ",
         "objects created with pairwise_interaction",
         call. = FALSE
    )
  }

  # pull out dynamics$hex for all and use to find unique objects
  dynamics <- get_unique_dynamics(dots)

  # create pairwise interaction matrix (binary matrix)
  hex_list <- sapply(dynamics, function(x) x$hex)
  structure <- matrix(0, nrow = length(hex_list), ncol = length(hex_list))
  for (i in seq_along(dots)) {
    target <- dots[[i]]$target$hex
    source <- dots[[i]]$source$hex
    row_id <- which(hex_list == target)
    col_id <- which(hex_list == source)
    structure[row_id, col_id] <- 1
  }

  # extract interaction effects for each pairwise interaction
  interaction_list <- lapply(dots, function(x) x$interaction)

  # and collate these into a list with one element for each species
  interaction <- vector("list", length = length(dynamics))
  interaction <- lapply(dynamics,
                        define_interaction,
                        dots = dots,
                        hex_list = hex_list,
                        interaction_list = interaction_list
  )

  # return
  as_multispecies(
    list(
      nspecies = length(dynamics),
      structure = structure,
      dynamics = dynamics,
      interaction = interaction
    )
  )

}

#' @name pairwise_interaction
#' @title Specify interactions between two species
#' @description Define population dynamics for multiple species from
#'   a set of single-species \code{\link{dynamics}} objects and
#'   defined pairwise interactions.
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
#' @details To be completed.
pairwise_interaction <- function(target, source, masks, funs) {

  # force evaluation to avoid NULL functions down the line
  force(masks)
  force(funs)

  # define function that specifies effects of source on target
  if (is.list(masks)) {
    interaction <- function(x, n, ...) {
      for (i in seq_along(masks)) {
        x <- do_mask(x, masks[[i]], funs[[i]], n, ...)
      }
      x
    }
  } else {
    interaction <- function(x, n, ...) {
      do_mask(x, masks, funs, n, ...)
    }
  }

  # return
  as_interaction(
    list(
      target = target,
      source = source,
      interaction = interaction
    )
  )

}

# internal function: extract unique dynamics objects from
#   list of interaction objects
get_unique_dynamics <- function(interactions) {
  # which hex values do we have?
  hex_target <- sapply(interactions, function(x) x$target$hex)
  hex_source <- sapply(interactions, function(x) x$source$hex)

  # create a vector of unique targets and sources
  hex_list <- unique(c(hex_target, hex_source))

  # which pops match these unique hexes?
  target_ids <- match(hex_list, hex_target)
  out <- lapply(
    interactions[target_ids[!is.na(target_ids)]],
    \(x) x$target
  )
  if (any(is.na(target_ids))) {
    source_ids <- match(hex_list[is.na(target_ids)], hex_source)
    out <- c(
      out,
      lapply(interactions[source_ids], \(x) x$source)
    )
  }

  # return all required dynamics objects
  out

}

# internal function: define interaction of each species with any other
#   species
define_interaction <- function(dyn, dots, hex_list, interaction_list) {
  # TODO: consider generalising to mean field interactions
  # aggregate_species <- user_defined_function(n, source) {
  #    n[source[1]][1:5] + n[source[2]][6:10]
  # }
  # source_sp <- match(sapply(dots[idx], function(x) x$source$hex), hex_list)
  # out <- function(x, n, ...) {
  #   for (j in seq_along(idx)) {
  #     x <- interaction_list[[idx[j]]](
  #       x, aggregate_species(n, source_sp), ...
  #     )
  #   }
  #   x
  # }

  # pull out relevant elements of fn_list
  idx <- which(sapply(dots, function(x) x$target$hex) == dyn$hex)

  # create function is species i is a target for any other species
  if (length(idx) > 0) {
    source_sp <- match(sapply(dots[idx], function(x) x$source$hex), hex_list)
    out <- function(x, n, ...) {
      for (j in seq_along(idx)) {
        x <- interaction_list[[idx[j]]](x, n[[source_sp[j]]], ...)
      }
      x
    }
  } else {
    out <- function(x, n, ...) {
      identity(x)
    }
  }

  # return
  out
}

# S3 method
#' @rdname multispecies
#'
#' @export
#'
#' @param x an object to pass to \code{is.multispecies}
# nolint start
is.multispecies <- function(x) {
  # nolint end
  inherits(x, "multispecies")
}

# S3 method
#' @export
# nolint start
print.multispecies <- function(x, ...) {
  # nolint end
  cat(paste0("Multispecies population dynamics object\n"))
}

# S3 method
#' @rdname multispecies
#' @export
# nolint start
is.interaction <- function(x) {
  # nolint end
  inherits(x, "interaction")
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
