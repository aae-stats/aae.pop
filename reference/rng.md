# Random number generators not available in existing R packages

Draw random numbers from unusual distributions, such as on the unit or
non-negative real line with known means and standard deviations.

## Usage

``` r
rmultiunit(
  n,
  mean,
  sd,
  Sigma = NULL,
  Omega = NULL,
  perfect_correlation = FALSE
)

rmultiunit_from_real(
  n,
  mean_real,
  sd_real = NULL,
  Sigma_chol = NULL,
  perfect_correlation = FALSE
)

runit_from_real(n, mean_real, sd_real)

runit(n, mean, sd)

unit_to_real(unit_mean, unit_sd)
```

## Arguments

- n:

  number of random draws to simulate. Each draw is a vector of values
  with length equal to `length(mean)` and `length(sd)` and the resulting
  output has `n` rows and `length(mean)` columns

- mean:

  vector of mean values on the unit scale

- sd:

  vector of positive standard deviations

- Sigma:

  optional covariance matrix with dimensions of `length(mean)` by
  `length(mean)` defining covariances between each pair of values in
  `mean`. Note that only the correlation structure is retained from
  `Sigma`, so that standard deviations are still required

- Omega:

  optional correlation matrix with dimensions of `length(mean)` by
  `length(mean)` defining correlations between each pair of values in
  `mean`

- perfect_correlation:

  `logical`, if `TRUE` and `Sigma` and `Omega` are `NULL`, then all
  values in each replicate (row) are perfectly correlated with known
  mean and standard deviation. If `FALSE`, then all values in each
  replicate are completely uncorrelated

- mean_real:

  vector of mean values converted to real-line equivalents

- sd_real:

  vector of standard deviations converted to real-line equivalents

- Sigma_chol:

  Cholesky decomposition of covariance matrix converted to real-line
  equivalent

- unit_mean:

  vector of mean values on the unit interval

- unit_sd:

  vector of standard deviations on the unit interval

## Value

a vector or matrix of random draws from the `r*unit` set of functions

## Details

The r\*unit family of functions support simulation of values on the unit
interval based on a known mean, sd, and correlation structure. `runit`
and `runit_from_real` are vectorised univariate functions, and
`rmultiunit` and `rmultiunit_from_real` are multivariate versions of
these same functions. `runit` and `rmultiunit` provide simulated values
on the unit line with specified means, standard deviations, and
correlation/covariance structure (in the case of `rmultiunit`).

The \*\_from_real versions of these functions are helpers that use
pre-transformed estimates of parameters on the real line, calculated
with `unit_to_real`. These functions are exported because
`unit_to_real`, called within `runit` and `rmultiunit`, is slow.
Separating this into a separate step allows less frequent calculations
of this transformation using function or dynamic versions of `args` in
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md).

`unit_to_real` converts means and standard deviations from their values
on the unit line to their equivalent values on the real line.

The use of the different versions of these functions is illustrated in
the Macquarie perch example on the package
\[website\](https://aae-stats.github.io/aae.pop/).

## Examples

``` r
# rmultiunit generates multivariate draws constrained to
#   the unit interval, with known mean, standard deviation,
#   and (optionally) covariance/correlation structure
rmultiunit(n = 10, mean = c(0.25, 0.5, 0.75), sd = c(0.1, 0.4, 0.25))
#>            [,1]       [,2]      [,3]
#>  [1,] 0.1575853 0.01761254 0.9445593
#>  [2,] 0.4276270 0.52709738 0.9603336
#>  [3,] 0.3186457 0.62405500 0.8450844
#>  [4,] 0.2722320 0.54031791 0.9504178
#>  [5,] 0.3991552 0.54087440 0.9689834
#>  [6,] 0.1838829 0.01500115 0.9838706
#>  [7,] 0.2273197 0.39881780 0.9967373
#>  [8,] 0.3448995 0.99999743 0.8872065
#>  [9,] 0.3162860 0.52977018 0.8180310
#> [10,] 0.2871536 0.96621961 0.8767448

# add in a correlation structure
omega_set <- cbind(
  c(1, 0.25, 0.01),
  c(0.25, 1, 0.5),
  c(0.01, 0.5, 1)
)
rmultiunit(
  n = 10,
  mean = c(0.25, 0.5, 0.75),
  sd = c(0.1, 0.4, 0.25),
  Omega = omega_set
)
#> Warning: NaNs produced
#> Warning: NaNs produced
#>            [,1]         [,2]       [,3]
#>  [1,] 0.1390801 8.584424e-01 0.88133202
#>  [2,] 0.2980836 9.850293e-01 0.92972847
#>  [3,] 0.2346260 9.989593e-01 0.99084645
#>  [4,] 0.1048095 2.155613e-04 0.09527154
#>  [5,] 0.3132903 9.935197e-01 0.99884698
#>  [6,] 0.2815907 1.366113e-06 0.13019641
#>  [7,] 0.4289172 9.999996e-01 0.97808043
#>  [8,] 0.1296939 6.571871e-01 0.53406620
#>  [9,] 0.2944634 5.203511e-04 0.36741513
#> [10,] 0.2610939 3.688165e-02 0.93892145
```
