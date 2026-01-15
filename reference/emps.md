# Calculate expected minimum population size (EMPS) for a [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md) object

Calculate expected minimum population size (EMPS) for a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object

## Usage

``` r
emps(sims, subset = NULL, times = NULL, fun = mean, ...)
```

## Arguments

- sims:

  an object returned from
  [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)

- subset:

  `integer` vector denoting the population classes to include in
  calculation of population abundance. Defaults to all classes

- times:

  `integer` vector specifying generations to include in calculation of
  extinction risk. Defaults to all simulated generations

- fun:

  `function` used to calculate average over all replicate trajectories.
  Defaults to `mean`. Alternatives might include `median` or `min`

- ...:

  additional arguments passed to `fun`

## Details

Expected minimum population size (EMPS) is the average minimum value of
all replicate trajectories. This value represents an expected lower
bound on population sizes over all generations, accounting for variation
among replicates. Abundances can be specified for all population classes
or for a subset of classes.

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
sims <- simulate(dyn, nsim = 1000)

# calculate expected minimum population size
emps(sims)
#> [1] 35.65322

# calculate expected minimum population size for 4 and 5 year
#   olds only
emps(sims, subset = 4:5)
#> [1] 1.320385

# calculate expected minimum population size but ignore first 10 years
emps(sims, times = 11:51)
#> [1] 35.99845

# calculate expected minimum population size based on median
emps(sims, fun = median)
#> [1] 35.6637
```
