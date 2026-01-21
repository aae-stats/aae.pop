# Calculate (quasi-)extinction risk at multiple thresholds for a [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md) object

Calculate (quasi-)extinction risk at multiple thresholds for a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object

## Usage

``` r
risk_curve(sims, threshold = NULL, subset = NULL, times = NULL, n = 100)
```

## Arguments

- sims:

  an object returned from
  [`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)

- threshold:

  `integer` or `numeric` vector denoting the set of threshold population
  sizes used to define the risk curve. Defaults to `n` evenly spaced
  values from 0 to the maximum observed abundance

- subset:

  `integer` vector denoting the population classes to include in
  calculation of population abundance. Defaults to all classes

- times:

  `integer` vector specifying generations to include in calculation of
  extinction risk. Defaults to all simulated generations

- n:

  `integer` specifying number of threshold values to use in default case
  when `threshold` is not specified. Defaults to 100

## Value

a named vector containing the threshold values (names) and the
probability the population will fall below these threshold values

## Details

Risk curves represent
[`pr_extinct`](https://aae-stats.github.io/aae.pop/reference/pr_extinct.md)
at multiple threshold population sizes simultaneously. This gives an
expression of risk of population declines below a range of values. Risk
curves are extracted from a
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
object as the proportion of replicate trajectories that fall below each
threshold value at any time step within a set period. Abundances can be
specified for all population classes or for a subset of classes.

The
[`get_cdf`](https://aae-stats.github.io/aae.pop/reference/get_cdf.md)
function is a much faster way to generate risk curves for almost all use
cases. The exception is when the `threshold` argument is used to specify
threshold values that are not evenly spaced.

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
sims <- simulate(dyn, nsim = 100)

# calculate risk curve
risk_curve(sims, n = 10)
#>   0  55 109 164 218 273 327 382 436 491 
#>   0   1   1   1   1   1   1   1   1   1 
```
