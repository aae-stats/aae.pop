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

# calculate risk curve
risk_curve(sims)
#>     0     6    12    18    25    31    37    43    49    55    61    67    74 
#> 0.000 0.000 0.000 0.000 0.020 0.199 0.564 0.888 0.982 1.000 1.000 1.000 1.000 
#>    80    86    92    98   104   110   116   123   129   135   141   147   153 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   159   165   172   178   184   190   196   202   208   215   221   227   233 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   239   245   251   257   264   270   276   282   288   294   300   306   313 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   319   325   331   337   343   349   355   362   368   374   380   386   392 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   398   405   411   417   423   429   435   441   447   454   460   466   472 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   478   484   490   496   503   509   515   521   527   533   539   546   552 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   558   564   570   576   582   588   595   601   607 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 

# calculate risk curve for multiple thresholds between 0 and 100
risk_curve(sims, threshold = seq(0, 100, by = 1))
#>     0     1     2     3     4     5     6     7     8     9    10    11    12 
#> 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
#>    13    14    15    16    17    18    19    20    21    22    23    24    25 
#> 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.002 0.006 0.008 0.009 0.015 0.025 
#>    26    27    28    29    30    31    32    33    34    35    36    37    38 
#> 0.042 0.065 0.091 0.126 0.171 0.228 0.276 0.325 0.400 0.473 0.516 0.578 0.646 
#>    39    40    41    42    43    44    45    46    47    48    49    50    51 
#> 0.703 0.762 0.821 0.855 0.891 0.924 0.944 0.958 0.972 0.979 0.981 0.987 0.993 
#>    52    53    54    55    56    57    58    59    60    61    62    63    64 
#> 0.996 0.999 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    65    66    67    68    69    70    71    72    73    74    75    76    77 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    78    79    80    81    82    83    84    85    86    87    88    89    90 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    91    92    93    94    95    96    97    98    99   100 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 

# calculate risk curve for 4 and 5 year olds only
risk_curve(sims, subset = 4:5)
#>     0     0     1     1     1     2     2     2     3     3     4     4     4 
#> 0.000 0.000 0.015 0.203 0.629 0.934 0.998 1.000 1.000 1.000 1.000 1.000 1.000 
#>     5     5     5     6     6     6     7     7     7     8     8     8     9 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>     9    10    10    10    11    11    11    12    12    12    13    13    13 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    14    14    14    15    15    16    16    16    17    17    17    18    18 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    18    19    19    19    20    20    21    21    21    22    22    22    23 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    23    23    24    24    24    25    25    25    26    26    27    27    27 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    28    28    28    29    29    29    30    30    30    31    31    31    32 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    32    33    33    33    34    34    34    35    35 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 

# calculate risk curve but ignore first 10 years
risk_curve(sims, times = 11:51)
#>     0     3     7    10    14    17    20    24    27    31    34    38    41 
#> 0.000 0.000 0.000 0.000 0.000 0.000 0.003 0.015 0.073 0.208 0.411 0.603 0.809 
#>    44    48    51    55    58    61    65    68    72    75    79    82    85 
#> 0.921 0.973 0.991 0.998 0.998 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>    89    92    96    99   102   106   109   113   116   120   123   126   130 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   133   137   140   143   147   150   154   157   161   164   167   171   174 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   178   181   184   188   191   195   198   202   205   208   212   215   219 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   222   225   229   232   236   239   243   246   249   253   256   260   263 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   266   270   273   277   280   284   287   290   294   297   301   304   307 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
#>   311   314   318   321   325   328   331   335   338 
#> 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
```
