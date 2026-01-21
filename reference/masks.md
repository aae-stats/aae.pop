# Isolate elements of population dynamics models

Helper functions to isolate particular components of a population
dynamics model, such as the reproduction terms, transition/growth terms,
or particular life stages from an abundance vector, such as pre- or
post-reproductive stages.

## Usage

``` r
reproduction(matrix, dims = NULL)

survival(matrix, dims = NULL)

transition(matrix, dims = NULL)

all_cells(matrix, dims = NULL)

all_classes(matrix, dims = NULL)

combine(...)
```

## Arguments

- matrix:

  a population dynamics matrix for which a particular mask is required.
  Only used to determine mask dimensions, so can be any matrix with
  appropriate dimensions

- dims:

  a numeric value or vector identifying subsets of cells to include in a
  given mask

- ...:

  a set of masks or masking functions to be combined into a single mask
  by one of the `combine` methods

## Value

`mask` object used to define the cells affected by a process included in
[`dynamics`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)

## Examples

``` r
# define a population
nclass <- 5
popmat <- matrix(0, nrow = nclass, ncol = nclass)
popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)

# pull out reproductive elements
reproduction(popmat)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] FALSE  TRUE  TRUE  TRUE  TRUE
#> [2,] FALSE FALSE FALSE FALSE FALSE
#> [3,] FALSE FALSE FALSE FALSE FALSE
#> [4,] FALSE FALSE FALSE FALSE FALSE
#> [5,] FALSE FALSE FALSE FALSE FALSE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# what if only 4 and 5 year olds reproduce?
reproduction(popmat, dims = 4:5)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] FALSE FALSE FALSE  TRUE  TRUE
#> [2,] FALSE FALSE FALSE FALSE FALSE
#> [3,] FALSE FALSE FALSE FALSE FALSE
#> [4,] FALSE FALSE FALSE FALSE FALSE
#> [5,] FALSE FALSE FALSE FALSE FALSE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# define survival elements
survival(popmat)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,]  TRUE FALSE FALSE FALSE FALSE
#> [2,] FALSE  TRUE FALSE FALSE FALSE
#> [3,] FALSE FALSE  TRUE FALSE FALSE
#> [4,] FALSE FALSE FALSE  TRUE FALSE
#> [5,] FALSE FALSE FALSE FALSE  TRUE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# what if 1 and 2 year olds always transition?
survival(popmat, dims = 3:5)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] FALSE FALSE FALSE FALSE FALSE
#> [2,] FALSE FALSE FALSE FALSE FALSE
#> [3,] FALSE FALSE  TRUE FALSE FALSE
#> [4,] FALSE FALSE FALSE  TRUE FALSE
#> [5,] FALSE FALSE FALSE FALSE  TRUE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# and transitions
transition(popmat)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] FALSE FALSE FALSE FALSE FALSE
#> [2,]  TRUE FALSE FALSE FALSE FALSE
#> [3,] FALSE  TRUE FALSE FALSE FALSE
#> [4,] FALSE FALSE  TRUE FALSE FALSE
#> [5,] FALSE FALSE FALSE  TRUE FALSE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# combine transitions and reproduction of 4 and 5 year olds
combine(reproduction(popmat, dims = 4:5), transition(popmat))
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] FALSE FALSE FALSE  TRUE  TRUE
#> [2,]  TRUE FALSE FALSE FALSE FALSE
#> [3,] FALSE  TRUE FALSE FALSE FALSE
#> [4,] FALSE FALSE  TRUE FALSE FALSE
#> [5,] FALSE FALSE FALSE  TRUE FALSE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# can also mask the population vector in this way
# pull out all classes
all_classes(popmat)
#>      [,1]
#> [1,] TRUE
#> [2,] TRUE
#> [3,] TRUE
#> [4,] TRUE
#> [5,] TRUE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 

# and just 3-5 year olds
all_classes(popmat, dims = 3:5)
#>       [,1]
#> [1,] FALSE
#> [2,] FALSE
#> [3,]  TRUE
#> [4,]  TRUE
#> [5,]  TRUE
#> attr(,"class")
#> [1] "mask"   "matrix" "array" 
```
