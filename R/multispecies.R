# combine dynamics in a multispecies object
multispecies <- function(interactions, ...) {

  # all ... should be dynamics objects

  # set outputs from input dynamics objects
  #  - decide how to handle conflicting settings for input pops
  list(
    ntime = defaults$ntime,
    nclass = nclass_list,
    nspecies = nspecies,
    matrix = matrix_list,
    stochasticity = stochasticity_list,
    density = density_list,
    rescale = rescale_list,
    is_modified = is_modified_list,
    is_stochastic = is_stochastic_list,
    is_density = is_density_list,
    is_rescaled = is_rescaled_list,
    is_multispecies = TRUE
  )

}
