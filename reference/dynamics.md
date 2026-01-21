# Create and update population dynamics objects

Define population dynamics from a matrix and additional objects that
determine covariate effects, density dependence, and forms of
stochasticity.

## Usage

``` r
dynamics(matrix, ...)

# S3 method for class 'dynamics'
update(object, ...)

is.dynamics(x)
```

## Arguments

- matrix:

  a matrix of vital rates specifying transitions between ages or stages.
  Specified in the format ntp1 = A A is the matrix and nt is the vector
  of abundances, so that values in a given column and row denote a
  transition from that column to that row

- ...:

  additional objects used to define population dynamics. Must be one or
  more of
  [`covariates`](https://aae-stats.github.io/aae.pop/reference/covariates.md),
  [`replicated_covariates`](https://aae-stats.github.io/aae.pop/reference/replicated_covariates.md),
  [`environmental_stochasticity`](https://aae-stats.github.io/aae.pop/reference/stochasticity.md),
  [`demographic_stochasticity`](https://aae-stats.github.io/aae.pop/reference/stochasticity.md),
  [`density_dependence`](https://aae-stats.github.io/aae.pop/reference/density_dependence.md),
  [`add_remove_pre`](https://aae-stats.github.io/aae.pop/reference/add_remove.md),
  or
  [`add_remove_post`](https://aae-stats.github.io/aae.pop/reference/add_remove.md).
  Note that
  [`density_dependence_n`](https://aae-stats.github.io/aae.pop/reference/density_dependence.md)
  is equivalent to
  [`add_remove_post`](https://aae-stats.github.io/aae.pop/reference/add_remove.md).

- object:

  a `dynamics` object

- x:

  an object to pass to `is.dynamics`

## Value

`dynamics` object containing a matrix population model and all
associated processes

## Details

A call to `dynamics` defines an object of class `dynamics`, which can be
used to simulate population trajectories with the
[`simulate`](https://aae-stats.github.io/aae.pop/reference/simulate.md)
function. The `plot` function is supported and will generate a general
life-cycle diagram based on the defined population dynamics.

A compiled `dynamics` object can be updated to change any of the
included processes with the `update` function.

## Examples

``` r
# define a population
nclass <- 5
popmat <- matrix(0, nrow = nclass, ncol = nclass)
popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)

# define a dynamics object
dyn <- dynamics(popmat)

# and plot this
if (rlang::is_installed("DiagrammeR")) {
  plot(dyn)
}

{"x":{"diagram":"digraph {\n\ngraph [layout = \"dot\",\n       outputorder = \"edgesfirst\",\n       bgcolor = \"white\",\n       rankdir = \"LR\"]\n\nnode [fontname = \"Helvetica\",\n      fontsize = \"10\",\n      shape = \"circle\",\n      fixedsize = \"true\",\n      width = \"0.5\",\n      style = \"filled\",\n      fillcolor = \"aliceblue\",\n      color = \"gray70\",\n      fontcolor = \"gray50\"]\n\nedge [fontname = \"Helvetica\",\n     fontsize = \"8\",\n     len = \"1.5\",\n     color = \"gray80\",\n     arrowsize = \"0.5\"]\n\n  \"1\" [label = \"Age 1\", fontcolor = \"#4F7942\", fontsize = \"12\", penwidth = \"2\", shape = \"circle\", color = \"#4F7942\", fillcolor = \"#4F794250\", width = \"0.9\", height = \"0.72\"] \n  \"2\" [label = \"Age 2\", fontcolor = \"#4F7942\", fontsize = \"12\", penwidth = \"2\", shape = \"circle\", color = \"#4F7942\", fillcolor = \"#4F794250\", width = \"0.9\", height = \"0.72\"] \n  \"3\" [label = \"Age 3\", fontcolor = \"#4F7942\", fontsize = \"12\", penwidth = \"2\", shape = \"circle\", color = \"#4F7942\", fillcolor = \"#4F794250\", width = \"0.9\", height = \"0.72\"] \n  \"4\" [label = \"Age 4\", fontcolor = \"#0E4D92\", fontsize = \"12\", penwidth = \"2\", shape = \"circle\", color = \"#0E4D92\", fillcolor = \"#0E4D9250\", width = \"0.9\", height = \"0.72\"] \n  \"5\" [label = \"Age 5\", fontcolor = \"#0E4D92\", fontsize = \"12\", penwidth = \"2\", shape = \"circle\", color = \"#0E4D92\", fillcolor = \"#0E4D9250\", width = \"0.9\", height = \"0.72\"] \n\"1\"->\"2\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"transition\", style = \"solid\"] \n\"2\"->\"3\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"transition\", style = \"solid\"] \n\"3\"->\"4\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"transition\", style = \"solid\"] \n\"4\"->\"1\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"reproduction\", style = \"dashed\"] \n\"4\"->\"5\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"transition\", style = \"solid\"] \n\"5\"->\"1\" [color = \"Gainsboro\", fontname = \"Avenir\", fontcolor = \"LightGray\", fontsize = \"14\", penwidth = \"4\", label = \"reproduction\", style = \"dashed\"] \n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
