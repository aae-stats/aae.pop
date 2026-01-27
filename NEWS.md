# aae.pop 0.2.0

## Features

- consistent handling of processes for all model types
- isolate extension packages and remove vignettes that depend on these

## Fixes

- address all checks for pre-CRAN submission

# aae.pop 0.1.2

## Features

- `add_remove_pre` and `add_remove_post` added to handle removals or additions
    to the population vector prior to or following the update step
- `plot` function for `dynamics` objects updated to allow custom labels

## Fixes

- Increase test coverage
- Update labelling for plots of dynamic objects to override all defaults
    when custom labels are supplied
- Update labelling for plots of dynamic objects to classify reproductive 
    first stages as reproduction, not survival
- Change `fecundity` to `reproduction` in `DiagrammeR` plots of dynamics objects
- Change `reproduction` mask to allow first stage to reproduce, with a default
    that excludes this first stage (`2:ncol(mat)`)
- Change `plot.dynamics` so that cycles in the first stage are labelled
    as reproduction by default. The argument `cycles_first` can be set
    to `"survival"` to revert to the original behaviour

# aae.pop 0.1.1

## Features

- `replicated_covariates` process class added to handle covariates
    applied separately to each replicate trajectory

# aae.pop 0.1.0

## Features

- `get_cdf`, `get_pdf`, `emps`, `exps` methods added to summarise
    simulated population trajectories: 

## Fixes

* replace Travis CI with GitHub Actions

## API changes

* removed deprecated `args.dyn` and `args.fn` arguments in `simulate`
* removed inference methods, now included in `aae.pop.inference` package
