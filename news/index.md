# Changelog

## aae.pop (development version)

## aae.pop 0.1.2

### Features

- New process class to handle removals or additions to the population
  vector prior to or following the update step: `add_remove_pre`,
  `add_remove_post`
- Allow custom labels on plotted dynamics objects with `labels` argument
- Update vignettes/website to include newer processes

### Fixes

- Increase test coverage
- Update labelling for plots of dynamic objects to override all defaults
  when custom labels are supplied
- Update labelling for plots of dynamic objects to classify reproductive
  first stages as reproduction, not survival
- Change `fecundity` to `reproduction` in `DiagrammeR` plots of dynamics
  objects
- Change `reproduction` mask to allow first stage to reproduce, with a
  default that excludes this first stage (`2:ncol(mat)`)
- Change `plot.dynamics` so that cycles in the first stage are labelled
  as reproduction by default. The argument `cycles_first` can be set to
  `"survival"` to revert to the original behaviour

## aae.pop 0.1.1

### Features

- New process class to handle covariates applied separately to each
  replicate trajectory: `replicated_covariates`

## aae.pop 0.1.0

### Features

- New methods to summarise simulated population trajectories: `get_cdf`,
  `get_pdf`, `emps`, `exps`

### Fixes

- replace Travis CI with GitHub Actions

### API changes

- removed deprecated `args.dyn` and `args.fn` arguments in `simulate`
- removed inference methods, now included in `aae.pop.inference` package
