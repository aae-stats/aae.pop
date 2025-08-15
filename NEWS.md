# aae.pop (development version)

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
