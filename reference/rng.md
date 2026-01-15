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
#>            [,1]         [,2]      [,3]
#>  [1,] 0.3258284 7.109279e-01 0.7253066
#>  [2,] 0.3147578 7.873457e-02 0.7816593
#>  [3,] 0.3254525 3.666655e-04 0.6590754
#>  [4,] 0.1735401 8.574929e-01 0.1742375
#>  [5,] 0.2198729 7.748146e-08 0.9541984
#>  [6,] 0.2834658 8.614973e-05 0.5670721
#>  [7,] 0.2843562 4.763717e-01 0.9420337
#>  [8,] 0.3666019 9.007615e-02 0.8181540
#>  [9,] 0.1918320 3.397002e-03 0.4181535
#> [10,] 0.2054899 9.985405e-01 0.8860253

if (FALSE) { # \dontrun{
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
} # }
```
