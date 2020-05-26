# expand dynamics into a metapop
metapopulation <- function(dispersal, ...) {

  # all ... should be dynamics objects

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
