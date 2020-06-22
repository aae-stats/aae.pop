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
         call. = FALSE)
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
  for (i in seq_along(dynamics)) {

    # pull out relevant elements of fn_list
    idx <- which(sapply(dots, function(x) x$target$hex) == dynamics[[i]]$hex)

    # create function is species i is a target for any other species
    if (length(idx) > 0) {
      source_sp <- match(sapply(dots[idx], function(x) x$source$hex), hex_list)
      interaction[[i]] <- function(x, n, ...) {
        for (j in seq_along(idx))
          x <- interaction_list[[idx[j]]](x, n[[source_sp[j]]], ...)
        x
      }
    } else {
      interaction[[i]] <- function(x, n, ...) {
        identity(x)
      }
    }

  }

  # pull out covariates, check them, and set timesteps if required
  covariate_idx <- sapply(dynamics, function(x) !is.null(x$covariates))
  include_covariates <- any(covariate_idx)
  ntime <- NULL
  if (include_covariates)
    ntime <- unique(sapply(dynamics[covariate_idx], function(x) x$covariates$ntime))

  # check covariates are either missing or included with same number of
  #   time steps for all species
  if (length(ntime) > 1) {
    stop("covariates must have same number of time steps for all ",
         "species (unless excluded) but ... includes species with ",
         clean_paste(ntime),
         " time steps",
         call. = FALSE)
  }

  # return
  as_multispecies(
    list(nspecies = length(dynamics),
         include_covariates = include_covariates,
         ntime = ntime,
         structure = structure,
         dynamics = dynamics,
         interaction = interaction)
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
#' @details To be completed.
pairwise_interaction <- function(target, source, masks, funs) {

  # define function that specifies effects of source on target
  if (is.list(masks)) {
    interaction <- function(x, n, ...) {
      for (i in seq_along(masks))
        x <- do_mask(x, masks[[i]], funs[[i]], n, ...)
      x
    }
  } else {
    interaction <- function(x, n, ...) {
      do_mask(x, masks, funs, n, ...)
    }
  }

  # return
  as_interaction(
    list(target = target,
         source = source,
         interaction = interaction)
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
    function(x) x$target
  )
  if (any(is.na(target_ids))) {
    source_ids <- match(hex_list[is.na(target_ids)], hex_source)
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
