# Calculate (quasi-)extinction risk for a [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md) object

Calculate (quasi-)extinction risk for a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object

## Usage

``` r
pr_extinct(sims, threshold = 0, subset = NULL, times = NULL)
```

## Arguments

- sims:

  an object returned from
  [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)

- threshold:

  `integer` or `numeric` denoting the threshold population size below
  which a population is considered functionally extinct. Defaults to `0`

- subset:

  `integer` vector denoting the population classes to include in
  calculation of population abundance. Defaults to all classes

- times:

  `integer` vector specifying generations to include in calculation of
  extinction risk. Defaults to all simulated generations

## Details

Quasi-extinction risk is the probability of decline below some specified
abundance threshold. This probability is extracted from a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object as the proportion of replicate trajectories that fall below this
threshold at any time step within a set period. Abundances can be
specified for all population classes or for a subset of classes.

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

# calculate quasi-extinction risk at a threshold population size
#   of 100 individuals
pr_extinct(sims, threshold = 100)
#> [1] 1

# repeat previous but focused on 4 and 5 year olds only
pr_extinct(sims, threshold = 100, subset = 4:5)
#> [1] 1

# repeat previous but ignore first 10 years
pr_extinct(sims, threshold = 100, times = 11:51)
#> [1] 1
```
