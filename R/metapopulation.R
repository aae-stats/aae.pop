#' @name metapopulation
#' @title Create a metapopulation dynamics object
#' @description Define population dynamics for multiple populations of
#'   a single species linked by dispersal (a metapopulation).
NULL

#' @rdname metapopulation
#'
#' @export
#'
#' @param structure binary or logical matrix denoting dispersal
#'   links between populations. Columns move to rows, so a \code{1} or
#'   \code{TRUE} in cell (a, b) denotes movement from population b
#'   to population a
#' @param dynamics a \code{\link{dynamics}} object or list of dynamics objects
#'   with one element for each population (each column/row of
#'   \code{structure}). If a single dynamics object is
#'   provided, it is recycled over all required populations
#' @param dispersal object created with
#'   \code{\link{dispersal}}. \code{dispersal} objects
#'   describe movements between populations and can include
#'   class-specific movements and density-dependent movements.
#'   \code{dispersal} objects must be a list with one element for each
#'   link in \code{structure}. These links are interpreted in
#'   column-major order, so that \code{dispersal} objects
#'   must be ordered by links in column 1, then column 2, and so on
#'
#' @details The \code{metapopulation} function connects multiple
#'   populations through known dispersal probabilities, handling
#'   standardisations of dispersal probabilities (if required) and
#'   updating masks and functions for all processes defined within each
#'   population. Further details on the definition of dispersal terms
#'   are provided in \code{\link{dispersal}}.
#'
#'   Covariates can be included in metapopulation models. The default
#'   behaviour is for all populations to share a single set of covariates,
#'   with covariate associations and masks defined separately for each
#'   population. A workaround to the assumption of shared covariates
#'   is included in the examples, below. Including covariates on dispersal
#'   probabilities requires covariate associations and masks defined
#'   on the combined metapopulation model. This approach is possible
#'   but currently untested.
#'
#' @returns \code{metapopulation} object containing a matrix metapopulation
#'   model; for use with \code{\link{simulate}}
#'
#' @examples
#' # define some populations, all with identical vital rates
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- lapply(
#'   1:3,
#'   function(i) dynamics(popmat)
#' )
#'
#' # define metapopulation structure with populations
#' #   1 and 3 dispersing into population 2
#' pop_structure <- matrix(0, nrow = 3, ncol = 3)
#' pop_structure[1, 2] <- 1
#' pop_structure[3, 2] <- 1
#'
#' # define dispersal between populations
#' dispersal_matrix <- matrix(0, nrow = nclass, ncol = nclass)
#' dispersal_matrix[survival(dispersal_matrix, dims = 20:25)] <- 0.2
#' pop_dispersal1 <- dispersal(dispersal_matrix, proportion = TRUE)
#' pop_dispersal2 <- dispersal(dispersal_matrix, proportion = FALSE)
#' pop_dispersal <- list(pop_dispersal1, pop_dispersal2)
#'
#' # create metapopulation object
#' metapop <- metapopulation(pop_structure, dyn, pop_dispersal)
#'
#' # simulate without covariates
#' sims <- simulate(metapop, nsim = 2)
#'
#' # simulate with shared covariates
#' # define a covariates function
#' covar_fn <- function(mat, x) {
#'   mat * (1 / (1 + exp(-0.5 * (x + 10))))
#' }
#'
#' # and some covariates
#' xsim <- matrix(rnorm(20), ncol = 1)
#'
#' # update the population dynamics objects with covariates
#' dyn <- lapply(
#'   dyn,
#'   update,
#'   covariates(masks = transition(dyn[[1]]$matrix), funs = covar_fn)
#' )
#'
#' # (re)create metapopulation object
#' metapop <- metapopulation(pop_structure, dyn, pop_dispersal)
#' sims <- simulate(
#'   metapop,
#'   nsim = 2,
#'   args = list(covariates = format_covariates(xsim))
#' )
#'
#' # simulate with separate covariates
#' #  (requires re-definition of covariate functions)
#' new_fn <- function(i) {
#'   force(i)
#'   function(mat, x) {
#'     mat * (1 / (1 + exp(-0.5 * (x[i] + 10))))
#'   }
#' }
#' new_fn <- lapply(
#'   1:3,
#'   new_fn
#' )
#'
#' # update the population dynamics objects with covariates
#' dyn <- lapply(
#'   dyn,
#'   update,
#'   covariates(masks = transition(dyn[[1]]$matrix), funs = covar_fn)
#' )
#'
#' # (re)create metapopulation object
#' metapop <- metapopulation(pop_structure, dyn, pop_dispersal)
#'
#' # and simulate with one column of predictors for each population
#' xsim <- matrix(rnorm(60), ncol = 3)
#' sims <- simulate(
#'   metapop,
#'   nsim = 2,
#'   args = list(covariates = format_covariates(xsim))
#' )
metapopulation <- function(structure, dynamics, dispersal) {
  # check the structure is ok
  structure <- check_structure(structure)

  # check dispersal is OK and convert to list if ndispersal = 1
  dispersal <- check_dispersal(dispersal, structure$ndispersal)

  # dynamics can be a list of length npop or a single matrix,
  #   expand if single matrix, error if multiple dynamics but not
  #   npop dynamics
  if (is.dynamics(dynamics)) {
    dynamics <- lapply(seq_len(structure$npop), function(x) dynamics)
  } else {
    if (length(dynamics) != structure$npop) {
      stop("dynamics must be a single dynamics object or a list of dynamics ",
           "objects with one element for each row/col of structure",
           call. = FALSE
      )
    }
  }

  # check list of dynamics objects is ok
  dyn_check <- check_dynamics(dynamics)

  # create block diagonal with dynamics matrices
  metapop_matrix <- block_diagonal(lapply(dynamics, function(x) x$matrix))

  # add dispersal to off-diagonal elements where structure == 1
  str_rows <- row(structure$structure)[structure$structure]
  str_cols <- col(structure$structure)[structure$structure]
  metapop_matrix <- add_dispersal(
    metapop_matrix, str_rows, str_cols, dispersal$dispersal, dyn_check$nclass
  )

  # define masks for each population
  mat_skeleton <- metapop_matrix
  if (is.list(mat_skeleton)) {
    mat_skeleton <- mat_skeleton[[1]]
  }
  mat_masks <- lapply(
    seq_len(structure$npop),
    function(i) metapop_idx(mat_skeleton, dyn_check$nclass, from = i, to = i)
  )
  pop_masks <- lapply(
    seq_len(structure$npop),
    function(i) ((i - 1) * dyn_check$nclass + 1):(i * dyn_check$nclass)
  )

  # define masks for each dispersal element
  dispersal_masks <- lapply(
    seq_len(structure$ndispersal),
    function(i) {
      metapop_idx(
        mat_skeleton,
        dyn_check$nclass,
        from = str_cols[i],
        to = str_rows[i]
      )
    }
  )

  # add in covariates if included in any populations
  covars <- NULL
  if (any(dyn_check$covars)) {
    # pull out functions for all populations
    covar_funs <- lapply(dynamics, function(x) x$covariates)

    # pull out non-NULL elements only
    covar_masks <- mat_masks[dyn_check$covars]
    covar_funs <- covar_funs[dyn_check$covars]

    # create full environmental stochasticity component
    covars <- covariates(covar_masks, covar_funs)
  }

  # add in covariates if included in any populations
  rep_covars <- NULL
  if (any(dyn_check$rep_covars)) {
    # pull out functions for all populations
    rep_covar_funs <- lapply(dynamics, function(x) x$replicated_covariates)

    # pull out non-NULL elements only
    rep_covar_masks <- mat_masks[dyn_check$rep_covars]
    rep_covar_funs <- rep_covar_funs[dyn_check$rep_covars]

    # create full environmental stochasticity component
    rep_covars <- covariates(rep_covar_masks, rep_covar_funs)
  }

  # add in environmental stochasticity if included in any populations
  envstoch <- NULL
  if (any(dyn_check$envstoch) || any(dispersal$stochasticity)) {
    # pull out functions for all populations
    env_funs <- lapply(dynamics, function(x) x$environmental_stochasticity)

    # pull out non-NULL elements only
    env_masks <- mat_masks[dyn_check$envstoch]
    env_funs <- env_funs[dyn_check$envstoch]

    # add in stochasticity for dispersal terms (if included)
    dispersal_stoch <- lapply(dispersal$dispersal, function(x) x$stochasticity)

    # add non-NULL elements to masks and funs
    missing <- sapply(dispersal_stoch, is.null)
    env_masks <- c(env_masks, dispersal_masks[!missing])
    env_funs <- c(env_funs, dispersal_stoch[!missing])

    # create full environmental stochasticity component
    envstoch <- environmental_stochasticity(env_masks, env_funs)
  }

  # expand demographic stochasticity if included
  demostoch <- NULL
  if (any(dyn_check$demostoch)) {
    # pull out masks and functions for all populations
    demo_funs <- lapply(dynamics, function(x) x$demographic_stochasticity)

    # pull out non-NULL elements only
    demo_masks <- pop_masks[dyn_check$demostoch]
    demo_funs <- demo_funs[dyn_check$demostoch]

    # create full demographic stochasticity component
    demostoch <- demographic_stochasticity(demo_masks, demo_funs)
  }

  # expand density dependence if included
  dens_depend <- NULL
  if (any(dyn_check$dens_depend) || any(dispersal$density)) {
    # pull out functions for all populations
    dens_funs <- lapply(dynamics, function(x) x$density_dependence)

    # pull out non-NULL elements only
    dens_masks <- mat_masks[dyn_check$dens_depend]
    nmasks <- pop_masks[dyn_check$dens_depend]
    dens_funs <- dens_funs[dyn_check$dens_depend]

    # add in stochasticity for dispersal terms (if included)
    dispersal_dens <- lapply(dispersal$dispersal, function(x) x$density)

    # add non-NULL elements to masks and funs
    missing <- sapply(dispersal_dens, is.null)
    dens_masks <- c(dens_masks, dispersal_masks[!missing])
    nmasks <- c(
      nmasks,
      lapply(
        seq_len(sum(!missing)),
        function(x) rep(TRUE, nrow(metapop_matrix))
      )
    )
    dens_funs <- c(dens_funs, dispersal_dens[!missing])

    # create full environmental stochasticity component
    dens_depend <- density_dependence(
      dens_masks,
      dens_funs,
      nmasks
    )
  }

  # expand rescale density dependence if included
  dens_depend_n <- NULL
  if (any(dyn_check$dens_depend_n)) {
    # pull out masks and functions for all populations
    dens_n_funs <- lapply(dynamics, function(x) x$density_dependence_n)

    # pull out non-NULL elements only
    dens_n_masks <- pop_masks[dyn_check$dens_depend_n]
    dens_n_funs <- dens_n_funs[dyn_check$dens_depend_n]

    # create full demographic stochasticity component
    dens_depend_n <- density_dependence_n(dens_n_masks, dens_n_funs)
  }

  # expand add_remove_pre if included
  add_rem_pre <- NULL
  if (any(dyn_check$add_rem_pre)) {
    # pull out masks and functions for all populations
    add_remove_pre_funs <- lapply(dynamics, function(x) x$add_remove_pre)

    # pull out non-NULL elements only
    add_remove_pre_masks <- pop_masks[dyn_check$add_rem_pre]
    add_remove_pre_funs <- add_remove_pre_funs[dyn_check$add_rem_pre]

    # create full demographic stochasticity component
    add_rem_pre <- add_remove_pre(add_remove_pre_masks, add_remove_pre_funs)
  }

  # expand add_remove_post if included
  add_rem_post <- NULL
  if (any(dyn_check$add_rem_post)) {
    # pull out masks and functions for all populations
    add_remove_post_funs <- lapply(dynamics, function(x) x$add_remove_post)

    # pull out non-NULL elements only
    add_remove_post_masks <- pop_masks[dyn_check$add_rem_post]
    add_remove_post_funs <- add_remove_post_funs[dyn_check$add_rem_post]

    # create full demographic stochasticity component
    add_rem_post <- add_remove_pre(add_remove_post_masks, add_remove_post_funs)
  }

  # collate metapop object with expanded dynamics
  metapop_dynamics <- list(
    nclass = dyn_check$nclass * structure$npop,
    npopulation = structure$npop,
    nspecies = 1,
    hex = hex_id(),
    matrix = metapop_matrix,
    covariates = covars,
    replicated_covariates = rep_covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = demostoch,
    density_dependence = dens_depend,
    density_dependence_n = dens_depend_n,
    add_remove_pre = add_rem_pre,
    add_remove_post = add_rem_post
  )

  # return metapopulation object
  as_metapopulation(metapop_dynamics)
}

# internal function: check the structure matrix is OK
check_structure <- function(x) {
  # is structure actually a matrix?
  if (!is.matrix(x)) {
    stop("structure must be a matrix", call. = FALSE)
  }

  # structure must be a square matrix
  if (nrow(x) != ncol(x)) {
    stop("structure must be a square matrix", call. = FALSE)
  }

  # convert binary matrix to logical
  if (all(x %in% c(0, 1))) {
    x <- x > 0
  }

  # is structure now logical?
  if (!all(x %in% c(TRUE, FALSE))) {
    stop(
      "structure must be binary (0/1) or logical (TRUE/FALSE)",
      call. = FALSE
    )
  }

  # remove diag if included
  if (any(diag(x))) {
    diag(x) <- FALSE
  }

  # how many populations are included?
  npop <- nrow(x)

  # how many dispersal elements exist?
  ndispersal <- sum(x[upper.tri(x)]) + sum(x[lower.tri(x)])

  # return list of key information
  list(
    npop = npop,
    ndispersal = ndispersal,
    structure = x
  )
}

# internal function: check dispersal object
check_dispersal <- function(x, n) {
  if (missing(x)) {
    stop(
      "dispersal must be provided to define a metapopulation",
      call. = FALSE
    )
  }

  # check dispersal is OK
  if (length(x) != n) {
    stop(
      "dispersal must have one element for each non-zero",
      " element of structure",
      call. = FALSE
    )
  }

  # and convert dispersal to list if needed
  if (n == 1 && class(x)[1] == "dispersal") {
    x <- list(x)
  }

  stochasticity <- check_processes(x, type = "stochasticity")
  density <- check_processes(x, type = "density")

  # return checked and formatted output
  list(
    dispersal = x, stochasticity = stochasticity, density = density
  )
}

# internal function: check list of dynamics object can be turned into
#   metapopulation object
check_dynamics <- function(dyn_list) {
  # do all elements have the same number of classes?
  classes <- sapply(dyn_list, function(x) x$nclass)
  if (length(unique(classes)) != 1) {
    stop(
      "all populations in dynamics must have the same number of classes",
      call. = FALSE
    )
  }

  # check covariates
  covars <- check_processes(dyn_list, type = "covariates")

  # check replicated covariates
  rep_covars <- check_processes(dyn_list, type = "replicated_covariates")

  # check environmental_stochasticity
  envstoch <- check_processes(dyn_list, type = "environmental_stochasticity")

  # check demographic_stochasticity
  demostoch <- check_processes(dyn_list, type = "demographic_stochasticity")

  # check density_dependence
  dens_depend <- check_processes(dyn_list, type = "density_dependence")

  # check density_dependence_n
  dens_depend_n <- check_processes(dyn_list, type = "density_dependence_n")

  # check add_remove_pre
  add_rem_pre <- check_processes(dyn_list, type = "add_remove_pre")

  # check add_remove_post
  add_rem_post <- check_processes(dyn_list, type = "add_remove_post")

  # return
  list(
    nclass = unique(classes),
    covars = covars,
    rep_covars = rep_covars,
    envstoch = envstoch,
    demostoch = demostoch,
    dens_depend = dens_depend,
    dens_depend_n = dens_depend_n,
    add_rem_pre = add_rem_pre,
    add_rem_post = add_rem_post
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

# internal function: check processes from multiple dynamics objects
check_processes <- function(x, type) {
  # pull out relevant process
  procs <- lapply(x, function(x) x[[type]])

  # return flag checking which pops include process
  !sapply(procs, is.null)
}

# internal function: add dispersal elements to metapopulation matrix
# TODO: check removing list mat option hasn't broken anything
#   Doesn't seem possible to end up with a list here given `block_diagnoal`
#   is used first
add_dispersal <- function(
    mat,
    str_rows,
    str_cols,
    dispersal,
    nclass
) {
  # and loop through all dispersals, updating metapop matrix one-by-one
  for (i in seq_along(dispersal)) {
    # loop if we have a list of matrices (i.e. if covariates are included)
    # REMOVED: cannot be a list
    # if (is.list(mat)) {
    #   # work out which cells we need to update
    #   #   (all matrices should have identical dims)
    #   idx <- metapop_idx(mat[[1]], nclass, from = str_cols[i], to = str_rows[i])
    #
    #   # and add in dispersal bits
    #   mat <- lapply(
    #     mat,
    #     do_mask,
    #     mask = idx,
    #     fun = function(x) dispersal[[i]]$kernel
    #   )
    # } else {
    # work out which cells we need to update
    idx <- metapop_idx(mat, nclass, from = str_cols[i], to = str_rows[i])

    # and add in dispersal bits
    mat <- do_mask(mat, mask = idx, function(x) dispersal[[i]]$kernel)
    # }
  }

  # rescale if needed
  rescale_needed <- sapply(dispersal, function(x) x$proportion)
  if (any(rescale_needed)) {
    col_sub <- str_cols[rescale_needed]
    row_sub <- str_rows[rescale_needed]
    # if (is.list(mat)) {
    #   mat <- lapply(mat, rescale_dispersal, nclass, col_sub, row_sub)
    # } else {
    mat <- rescale_dispersal(mat, nclass, col_sub, row_sub)
    # }
  }

  # and check survival doesn't exceed 1
  for (i in seq_along(dispersal)) {
    # if (is.list(mat)) {
    #   survival_issue <- lapply(
    #     seq_along(mat),
    #     function(j) {
    #       check_survival(
    #         mat[[j]], nclass, str_cols[i], idx,
    #         timestep = j
    #       )
    #     }
    #   )
    #   issue <- sapply(survival_issue, function(x) x$issue)
    #   classes <- unlist(lapply(survival_issue, function(x) x$classes))
    #   if (any(issue)) {
    #     message(
    #       "Survival (including dispersal) exceeds 1 for classes ",
    #       clean_paste(unique(classes)),
    #       " in timesteps ", clean_paste(which(issue)),
    #       " for population ", i
    #     )
    #   }
    # } else {
    survival_issue <- check_survival(mat, nclass, str_cols[i], idx)
    if (survival_issue$issue) {
      message(
        "Survival (including dispersal) exceeds 1 for classes ",
        clean_paste(survival_issue$classes)
      )
      # }
    }
  }

  # return
  mat
}

# internal function: rescale dispersal kernel and source matrices
rescale_dispersal <- function(mat, nclass, cols, rows) {
  # work out how many dispersing from each population
  n_disperse <- table(cols)
  dispersers <- as.numeric(names(n_disperse))

  # loop over each source population with active dispersers
  for (i in seq_along(dispersers)) {
    # pull out the source population
    idx <- metapop_idx(mat, nclass, dispersers[i], dispersers[i])
    source <- matrix(mat[idx], ncol = nclass)

    # pull out the dispersing proportions, setting non-target rows to zero
    col_subset <-
      ((dispersers[i] - 1) * nclass + 1):(dispersers[i] * nclass)
    idy <- col(mat) %in% col_subset
    rows_sub <- rows[cols == dispersers[i]]
    target_rows <- NULL
    for (i in seq_along(rows_sub)) {
      target_rows <- c(
        target_rows,
        ((rows_sub[i] - 1) * nclass + 1):(rows_sub[i] * nclass)
      )
    }
    idz <- !(row(mat) %in% target_rows) & col(mat) %in% col_subset
    off_target <- mat[idz]
    kernels <- mat
    kernels[idz] <- 0
    kernels <- matrix(kernels[idy], ncol = nclass)

    # work out proportion available, excluding fecundity
    idr <- options()$aae.pop_reproduction_mask(source)
    reprod <- source[idr]
    source[idr] <- 0
    available <- apply(source, 2, sum)

    # work out proportion leaving
    leave <- apply(kernels, 2, sum)

    # work out proportion staying in source
    remain <- 1 - leave

    # update each accordingly
    kernels <- sweep(kernels, 2, available, "*")
    source <- sweep(source, 2, remain, "*")

    # add fecundity back in
    source[idr] <- reprod

    # and update matrix, kernels first because they have zeros
    #   in the source population elements
    mat[idy] <- kernels
    mat[idz] <- off_target
    mat[idx] <- source
  }

  mat
}

# internal function: define masks for metapopulations based on pop IDs
metapop_idx <- function(mat, nclass, from, to) {
  # work out which rows correspond to "to"
  row_subset <-
    ((to - 1) * nclass + 1):(to * nclass)

  # work out which cols correspond to "from"
  col_subset <-
    ((from - 1) * nclass + 1):(from * nclass)

  # return cell subset
  row(mat) %in% row_subset & col(mat) %in% col_subset
}

# internal function: check implied survival with dispersal
check_survival <- function(mat, nclass, col, idx, timestep = NULL) {
  # pull out the from population, add dispersal, and check proportion surviving
  idy <- metapop_idx(mat, nclass, from = col, to = col)
  total <- matrix(mat[idy] + mat[idx], ncol = nclass)
  total[options()$aae.pop_reproduction_mask(total)] <- 0
  total_survival <- apply(total, 2, sum)

  # return TRUE to signal an issue with classes recorded
  list(
    issue = any(total_survival > 1),
    classes = which(total_survival > 1)
  )
}

# S3 method
#' @rdname metapopulation
#'
#' @export
#'
#' @param x an object to pass to \code{is.metapopulation}
# nolint start
is.metapopulation <- function(x) {
  # nolint end
  inherits(x, "metapopulation")
}

# S3 method
#' @export
# nolint start
print.metapopulation <- function(x, ...) {
  # nolint end
  cat(paste0(
    "Metapopulation dynamics object with ",
    x$npopulation,
    " populations\n"
  ))
}

# internal function: set metapopulation class
as_metapopulation <- function(x) {
  # it should just be a list coming in but we want it to be a
  #   dynamics object. This will check it's a list and add
  #   the dynamics class.
  x <- as_dynamics(x)

  # return
  as_class(x, name = "metapopulation", type = "dynamics")
}
