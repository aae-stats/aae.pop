# Functions for a single time-step update

Define how population abundances are updated from one time step to the
next. Functions can take any form but will only be vectorised across
replicates in limited situations.

## Usage

``` r
update_crossprod(pop, mat)

update_binomial_leslie(pop, mat)

update_multinomial(pop, mat)
```

## Arguments

- pop:

  current state of the population

- mat:

  matrix of vital rates used to update population state

## Value

a matrix containing population abundances in each stage of a matrix
population model. Contains one row for each replicate population
trajectory and one column for each population stage

## Details

Updaters can be changed through the `options` argument to `simulate` and
can also be changed globally for an R session by changing the global
option with, e.g., `options(aae.pop_update = update_binomial_leslie)`

`update_crossprod` updates abundances with a direct matrix
multiplication that does not include any form of demographic
stochasticity. This is the fastest update option and will vectorise
across replicates if the population matrix is not expanded by
[`environmental_stochasticity`](https://aae-stats.github.io/aae.pop/reference/stochasticity.md)
or
[`density_dependence`](https://aae-stats.github.io/aae.pop/reference/density_dependence.md).

`update_binomial_leslie` updates abundances with a direct RNG draw that
combines update with demographic stochasticity, assuming a Leslie
matrix.

`update_multinomial` updates abundances with a direct RNG draw that
combines update with demographic stochasticity, allowing for general
matrix forms (slower than `update_binomial_leslie`).

## Examples

``` r
# define a basic population
nstage <- 5
popmat <- matrix(0, nrow = nstage, ncol = nstage)
popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)

# define a dynamics object
dyn <- dynamics(popmat)

# simulate with the default updater
sims <- simulate(dyn)

# simulate with a multinomial updater
sims <- simulate(dyn, options = list(update = update_multinomial))
```
