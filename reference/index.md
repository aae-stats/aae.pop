# Package index

## Defining population dynamics models

- [`dynamics()`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)
  [`update(`*`<dynamics>`*`)`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)
  [`is.dynamics()`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)
  : Create and update population dynamics objects

## Simulating population dynamics

- [`simulate(`*`<dynamics>`*`)`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
  [`is.simulation()`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
  [`is.simulation_list()`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
  : Simulate single or multispecies population dynamics in R

## Summarising population dynamics

- [`pr_extinct()`](https://aae-stats.github.io/aae.pop/reference/pr_extinct.md)
  :

  Calculate (quasi-)extinction risk for a
  [`simulate`](https://rdrr.io/r/stats/simulate.html) object

- [`emps()`](https://aae-stats.github.io/aae.pop/reference/emps.md) :

  Calculate expected minimum population size (EMPS) for a
  [`simulate`](https://rdrr.io/r/stats/simulate.html) object

- [`risk_curve()`](https://aae-stats.github.io/aae.pop/reference/risk_curve.md)
  :

  Calculate (quasi-)extinction risk at multiple thresholds for a
  [`simulate`](https://rdrr.io/r/stats/simulate.html) object

- [`exps()`](https://aae-stats.github.io/aae.pop/reference/exps.md) :

  Calculate expected population size for a
  [`simulate`](https://rdrr.io/r/stats/simulate.html) object based on
  generic functions (ExPS)

- [`get_pdf()`](https://aae-stats.github.io/aae.pop/reference/get_pdf.md)
  :

  Calculate the probability density of a summary statistic across all
  iterations of a [`simulate`](https://rdrr.io/r/stats/simulate.html)
  object

- [`get_cdf()`](https://aae-stats.github.io/aae.pop/reference/get_cdf.md)
  :

  Calculate the cumulative distribution function of a summary statistic
  across all iterations of a
  [`simulate`](https://rdrr.io/r/stats/simulate.html) object

## Adding density dependence

- [`density_dependence()`](https://aae-stats.github.io/aae.pop/reference/density_dependence.md)
  [`density_dependence_n()`](https://aae-stats.github.io/aae.pop/reference/density_dependence.md)
  : Specify density dependence in models of population dynamics
- [`beverton_holt()`](https://aae-stats.github.io/aae.pop/reference/density_functions.md)
  [`ricker()`](https://aae-stats.github.io/aae.pop/reference/density_functions.md)
  : Common forms of density dependence

## Altering the population vector

- [`add_remove_pre()`](https://aae-stats.github.io/aae.pop/reference/add_remove.md)
  [`add_remove_post()`](https://aae-stats.github.io/aae.pop/reference/add_remove.md)
  : Specify additions or removals in models of population dynamics

## Adding stochasticity

- [`environmental_stochasticity()`](https://aae-stats.github.io/aae.pop/reference/stochasticity.md)
  [`demographic_stochasticity()`](https://aae-stats.github.io/aae.pop/reference/stochasticity.md)
  : Specify environmental and demographic stochasticity in models of
  population dynamics

## Adding covariates

- [`covariates()`](https://aae-stats.github.io/aae.pop/reference/covariates.md)
  [`format_covariates()`](https://aae-stats.github.io/aae.pop/reference/covariates.md)
  : Specify covariate dependence in models of population dynamics
- [`replicated_covariates()`](https://aae-stats.github.io/aae.pop/reference/replicated_covariates.md)
  : Specify replicate-specific covariate dependence in models of
  population dynamics

## Subsetting population matrices

- [`reproduction()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  [`survival()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  [`transition()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  [`all_cells()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  [`all_classes()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  [`combine()`](https://aae-stats.github.io/aae.pop/reference/masks.md)
  : Isolate elements of population dynamics models

## Updating population abundances

- [`update_crossprod()`](https://aae-stats.github.io/aae.pop/reference/updaters.md)
  [`update_binomial_leslie()`](https://aae-stats.github.io/aae.pop/reference/updaters.md)
  [`update_multinomial()`](https://aae-stats.github.io/aae.pop/reference/updaters.md)
  : Functions for a single time-step update

## Adding metapopulation structure

- [`metapopulation()`](https://aae-stats.github.io/aae.pop/reference/metapopulation.md)
  [`is.metapopulation()`](https://aae-stats.github.io/aae.pop/reference/metapopulation.md)
  : Create a metapopulation dynamics object
- [`dispersal()`](https://aae-stats.github.io/aae.pop/reference/dispersal.md)
  : Specify dispersal between populations in a metapopulation model

## Modelling multiple species

- [`multispecies()`](https://aae-stats.github.io/aae.pop/reference/multispecies.md)
  [`is.multispecies()`](https://aae-stats.github.io/aae.pop/reference/multispecies.md)
  [`is.interaction()`](https://aae-stats.github.io/aae.pop/reference/multispecies.md)
  : Create a population dynamics object with multiple species
- [`pairwise_interaction()`](https://aae-stats.github.io/aae.pop/reference/pairwise_interaction.md)
  : Specify interactions between two species

## Helpers

- [`rmultiunit()`](https://aae-stats.github.io/aae.pop/reference/rng.md)
  [`rmultiunit_from_real()`](https://aae-stats.github.io/aae.pop/reference/rng.md)
  [`runit_from_real()`](https://aae-stats.github.io/aae.pop/reference/rng.md)
  [`runit()`](https://aae-stats.github.io/aae.pop/reference/rng.md)
  [`unit_to_real()`](https://aae-stats.github.io/aae.pop/reference/rng.md)
  : Random number generators not available in existing R packages
