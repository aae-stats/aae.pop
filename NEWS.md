# aae.pop (development version)

## Features

- New process class to handle removals or additions to the population vector
    prior to or following the update step: `add_remove_pre`, `add_remove_post`

## Fixes

- Change `fecundity` to `reproduction` in `DiagrammeR` plots of dynamics objects

# aae.pop 0.1.1

## Features

- New process class to handle covariates applied separately to each replicate
    trajectory: `replicated_covariates`

# aae.pop 0.1.0

## Features

- New methods to summarise simulated population trajectories: `get_cdf`, `get_pdf`, `emps`, `exps`

## Fixes

* replace Travis CI with GitHub Actions

## API changes

* removed deprecated `args.dyn` and `args.fn` arguments in `simulate`
* removed inference methods, now included in `aae.pop.inference` package
