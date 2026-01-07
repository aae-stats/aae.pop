# Specify interactions between two species

Define population dynamics for multiple species from a set of
single-species
[`dynamics`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)
objects and defined pairwise interactions.

## Usage

``` r
pairwise_interaction(target, source, masks, funs)
```

## Arguments

- target:

  population whose vital rates are affected by the pairwise interaction

- source:

  population whose abundances affect the vital rates of `target`

- masks:

  masks defining which vital rates are influenced by each function

- funs:

  functions that take vital rates and abundances of the `source`
  population as inputs and return scaled vital rates

## Details

To be completed.
