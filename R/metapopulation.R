# expand dynamics into a metapop
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
  list(
    ntime = defaults$ntime,
    nclass = defaults$nclass,
    nspecies = defaults$nspecies,
    matrix = matrix,
    stochasticity = stochasticity,
    density = density,
    rescale = rescale,
    is_modified = defaults$is_modified,
    is_stochastic = defaults$is_stochastic,
    is_density = defaults$is_density,
    is_rescaled = defaults$is_rescaled,
    is_multispecies = FALSE
  )

}
