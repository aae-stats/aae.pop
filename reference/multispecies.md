# Create a population dynamics object with multiple species

Define population dynamics for multiple species from a set of
single-species
[`dynamics`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)
objects and defined pairwise interactions.

## Usage

``` r
multispecies(...)

is.multispecies(x)

is.interaction(x)
```

## Arguments

- ...:

  `pairwise_interaction` objects defining a set of pairwise interactions
  between species

- x:

  an object to pass to `is.multispecies`
