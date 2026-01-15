# Calculate expected population size for a [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md) object based on generic functions (ExPS)

Calculate expected population size for a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object based on generic functions (ExPS)

## Usage

``` r
exps(
  sims,
  subset = NULL,
  times = NULL,
  fun_within = mean,
  fun_among = mean,
  ...
)
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

- fun_within:

  `function` used to summarise a single trajectory. Must return a single
  value. Defaults to `mean`

- fun_among:

  `function` used to summarise over all replicate trajectories. Defaults
  to `mean`. Alternatives might include `median` or `min`

- ...:

  additional arguments passed to `fun_within` and `fun_among`. If these
  conflict, a wrapper function could be used to define expected
  arguments for each function

## Details

Expected population size (ExPS) is a highly flexible generalisation of
[`emps`](https://aae-stats.github.io/aae.pop/reference/emps.md) and
represents a two-level summary that first summarises individual
population trajectories and then summarises these values over all
replicates. Abundances can be specified for all population classes or
for a subset of classes.

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

# calculate expected population size
exps(sims)
#> [1] 96.072

# calculate expected population size for 4 and 5 year
#   olds only
exps(sims, subset = 4:5)
#> [1] 4.860821

# calculate expected population size but ignore first 10 years
exps(sims, times = 11:51)
#> [1] 79.90515

# calculate expected population size based on median
exps(sims, fun_among = median)
#> [1] 95.3747

# calculate expected maximum population size based on median
exps(sims, fun_within = max, fun_among = median)
#> [1] 322.7

# calculate exps with conflicting quantile functions, handling
#   conflicting arguments with wrapper functions
quant1 <- function(x, p1, ...) {
  quantile(x, prob = p1)
}
quant2 <- function(x, p2, ...) {
  quantile(x, prob = p2)
}
exps(
  sims,
  fun_within = quant1, fun_among = quant2, p1 = 0.25, p2 = 0.75
)
#>     75% 
#> 60.7253 
```
