#' @name metapopulation
#' @title Create a population dynamics object with multiple populations
#' @description Define population dynamics for multiple populations of
#'   a single species linked by dispersal.
NULL

#' @rdname metapopulation
#'
#' @export
#'
#' @param structure djkfd
#' @param dynamics dkjfd
#' @param dispersal dkfjd
#' @param \dots djfkd
#'
#' @details djkfdj
#'
#' @examples
#' # do stuff
metapopulation <- function(structure, dynamics, dispersal, ...) {

  # check the structure is ok
  structure <- check_structure(structure)

  # check dispersal is provided
  if (missing(dispersal))
    stop("dispersal must be provided to define a metapopulation", call. = FALSE)

  # check dispersal is OK
  if (length(dispersal) != structure$ndispersal)
    stop("dispersal must have one element for each non-zero element of structure", call. = FALSE)

  # and convert dispersal to list if needed
  if (structure$ndispersal == 1 & class(dispersal)[1] == "dispersal")
    dispersal <- list(dispersal)

  # dynamics can be a list of length npop or a single matrix,
  #   expand if single matrix
  if (class(dynamics)[1] == "dynamics") {
    dynamics <- lapply(seq_len(structure$npop), function(x) dynamics)
  }

  # check list of dynamics objects is ok
  dyn_check <- check_dynamics(dynamics)

  # expand environmental stochasticity if included
  covars <- NULL
  if (dyn_check$covars) {

    # expand, account for some without if that's the case
    covars <- lapply(dynamics, function(x) x$covariates)

  }

  # expand environmental stochasticity if included
  envstoch <- NULL
  if (dyn_check$envstoch) {

    # expand, account for some without if that's the case

    # maybe just need to make up new masks (one for each pop)
    #  and then set masks = masks and funs = lapply(dynamics, function(x) x$environmental_stochasticity)

    # include dispersal$stochasticity

  }

  # expand demographic stochasticity if included
  demostoch <- NULL
  if (dyn_check$demostoch) {

    # expand, account for some without if that's the case

  }

  # expand density dependence if included
  dens_depend <- NULL
  if (dyn_check$dens_depend) {

    # expand, account for some without if that's the case

    # include dispersal$density

  }

  # expand rescale density dependence if included
  dens_depend_n <- NULL
  if (dyn_check$dens_depend_n) {

    # expand, account for some without if that's the case

  }

  # create block diagonal with dynamics matrices
  if (dyn_check$covars) {
    metapop_matrix <- lapply(
      seq_len(dynamics$ntime),
      function(i) block_diagonal(lapply(dynamics, function(x) x$matrix[[i]]))
    )
  } else {
    metapop_matrix <- block_diagonal(lapply(dynamics, function(x) x$matrix))
  }

  # add dispersal to off-diagonal elements where structure == 1
  metapop_matrix <- add_dispersal(
    metapop_matrix, structure, dispersal, dyn_check$nstage
  )

  # collate metapop object with expanded dynamics
  metapop_dynamics <- list(
    ntime = dyn_check$ntime,
    nclass = nrow(metapop_matrix),
    npopulation = structure$npop,
    nspecies = 1,
    hex = hex_id(),
    base_matrix = lapply(dynamics, function(x) x$base_matrix),
    matrix = metapop_matrix,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = demostoch,
    density_dependence = dens_depend,
    density_dependence_n = dens_depend_n
  )

  # return metapopulation object
  as_metapopulation(metapop_dynamics)

}

# internal function: check the structure matrix is OK
check_structure <- function(x) {

  # is structure actually a matrix?
  if (!is.matrix(x))
    stop("structure must be a matrix", call. = FALSE)

  # convert binary matrix to logical
  if (all(x %in% c(0, 1)))
    x <- x > 0

  # is structure now logical?
  if (!all(x %in% c(TRUE, FALSE)))
      stop("structure must be binary (0/1) or logical (TRUE/FALSE)", call. = FALSE)

  # remove diag if included
  if (any(diag(x)))
    diag(x) <- FALSE

  # how many populations are included?
  npop <- nrow(x)

  # how many dispersal elements exist?
  ndispersal <- sum(x[upper.tri(x)]) + sum(x[lower.tri(x)])

  # return list of key information
  list(npop = npop,
       ndispersal = ndispersal,
       structure = x)

}

# internal function: check list of dynamics object can be turned into
#   metapopulation object
check_dynamics <- function(dyn_list) {

  # do all elements have the same number of stages?
  stages <- sapply(dyn_list, function(x) x$nclass)
  if (length(unique(stages)) != 1)
    stop("all populations in dynamics must have the same number of stages", call. = FALSE)

  # pull out covariates
  covars <- lapply(dyn_list, function(x) x$covariates)

  # are there any at all?
  if (all(sapply(covars, is.null))) {
    covars <- FALSE
    ntime <- 1
  } else {

    # if so, check they all have the same dimensions
    ncovar <- sapply(covars, function(x) x$ntime)
    if (length(unique(ncovar)) != 1)
      stop("all populations in dynamics must have the same number of time steps", call. = FALSE)

    # then if some are NULL and some are not, expand NULL matrices to lists
    if (any(sapply(covars, is.null))) {
      missing_covars <- sapply(covars, is.null)
      for (i in seq_along(dyn_list[missing_covars])) {
        dyn_list[[i]]$matrix <- lapply(seq_len(unique(covars)), function(x) dyn_list[[i]]$matrix)
        dyn_list[[i]]$ntime <- unique(ncovar)
      }
    }

    # set flag so we know we're dealing with lists of matrices
    covars <- TRUE
    ntime <- unique(ncovar)

  }

  # check environmental_stochasticity
  envstoch <- FALSE

  # check demographic_stochasticity
  demostoch <- FALSE

  # check density_dependence
  dens_depend <- FALSE

  # check density_dependence_n
  dens_depend_n <- FALSE

  # return
  list(
    ntime = ntime,
    nstage = unique(stages),
    covars = covars,
    envstoch = envstoch,
    demostoch = demostoch,
    dens_depend = dens_depend,
    dens_depend_n = dens_depend_n
  )

}

# internal function: create a block diagonal matrix from list of matrices
block_diagonal <- function(mats) {

  # work out dims and initialise an empty matrix
  nrow_single <- nrow(mats[[1]])
  nrows <- nrow_single * length(mats)
  mat <- matrix(0, nrow = nrows, ncol = nrows)

  # fill based on block locations
  for (i in seq_along(mats)) {
    subset <- ((i - 1) * nrow_single + 1):(i * nrow_single)
    idx <- row(mat) %in% subset & col(mat) %in% subset
    mat[idx] <- mats[[i]]
  }

  # return
  mat

}

# internal function: add dispersal elements to metapopulation matrix
add_dispersal <- function(mat, structure, dispersal, nstage) {

  # work out the populations corresponding to each dispersal
  str_rows <- row(structure$structure)[structure$structure]
  str_cols <- col(structure$structure)[structure$structure]

  # and loop through all dispersals, updating metapop matrix one-by-one
  for (i in seq_len(structure$ndispersal)) {

    # work out which cells we need to update
    idx <- metapop_idx(mat, nstage, from = str_cols[i], to = str_rows[i])

    # loop if we have a list of matrices (i.e. if covariates are included)
    if (is.list(mat)) {
      mat <- lapply(mat, do_mask, mask = idx, fun = function(x) dispersal[[i]]$kernel)
    } else {
      mat <- do_mask(mat, mask = idx, function(x) dispersal[[i]]$kernel)
    }

  }

  # return
  mat

}

# internal function: define masks for metapopulations based on pop IDs
metapop_idx <- function(mat, nstage, from, to) {

  # work out which rows correspond to "to"
  row_subset <-
    ((to - 1) * nstage + 1):(to * nstage)

  # work out which cols correspond to "from"
  col_subset <-
    ((from - 1) * nstage + 1):(from * nstage)

  # return cell subset
  row(mat) %in% row_subset & col(mat) %in% col_subset

}

# internal function: set metapopulation class
as_metapopulation <- function(x) {
  as_class(x, name = "metapopulation", type = "list")
}
