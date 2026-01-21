# Specify density dependence in models of population dynamics

Specify density dependence in vital rates (`density_dependence`) and in
total abundances (`density_dependence_n`).

## Usage

``` r
density_dependence(masks, funs, nmask = NULL)

density_dependence_n(masks, funs)
```

## Arguments

- masks:

  a logical matrix or vector (or list of these) defining cells affected
  by `funs`. See Details and
  [`masks`](https://aae-stats.github.io/aae.pop/reference/masks.md)

- funs:

  a function or list of functions with one element for each element of
  `masks`. See Details

- nmask:

  logical vector or list of vectors defining elements of the population
  vector affected by each mask-function pair. Intended primarily for
  internal use when scaling up processes in
  [`metapopulation`](https://aae-stats.github.io/aae.pop/reference/metapopulation.md)

## Value

`density_dependence` object specifying covariate effects on a matrix
population model; for use with
[`dynamics`](https://aae-stats.github.io/aae.pop/reference/dynamics.md)

## Details

`density_dependence` specifies standard density dependence on vital
rates, such as scramble or contest competition or allee effects.

`density_dependence_n` is an alternative parameterisation of density
dependence that acts directly on population abundances. Note that
`density_dependence_n` has been superseded by
[`add_remove_post`](https://aae-stats.github.io/aae.pop/reference/add_remove.md).

Masks must be of the same dimension as the population dynamics matrix
and specify cells influenced by density dependence according to `funs`.
In the case of `density_dependence_n`, `masks` are logical vectors with
one element for each class. Additional details on masks are provided in
[`masks`](https://aae-stats.github.io/aae.pop/reference/masks.md).

If using `density_depenence`, functions must take at least two
arguments, a matrix `x` and a vector `n`, which represent the population
dynamics matrix and the population abundances. Functions must return a
matrix with the same dimensions as `x`, modified to reflect the effects
of current abundances (`n`) on vital rates.

In the case of `density_dependence_n`, `funs` takes only one argument,
the population abundances `n` following all other updates in a given
iteration/generation. This allows rescaling of population abundances
based on total abundance or through more complicated functions that
depend on external arguments (e.g., mass mortality events or
harvesting).

Additional arguments to functions are supported and can be passed to
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
with the `args`, `args.dyn`, or `args.fn` arguments.

## Examples

``` r
# define a population matrix (columns move to rows)
nclass <- 5
popmat <- matrix(0, nrow = nclass, ncol = nclass)
popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)

# define a dynamics object
dyn <- dynamics(popmat)

# add some density dependence
dd <- density_dependence(
  masks = reproduction(popmat, dims = 4:5),
  funs = ricker(1000)
)

# update the dynamics object
dyn <- update(dyn, dd)

# simulate trajectories
sims <- simulate(dyn, nsim = 100, options = list(ntime = 50))

# and plot
plot(sims)
```
