## aae.pop: population dynamics models in R

aae.pop is a plug-and-play tool to simulate matrix population models.
It’s designed to be fast, flexible, and easily adapted to different
model structures, including [complex demographic
processes](https://aae-stats.github.io/aae.pop/articles/including_processes.md),
[metapopulations](https://aae-stats.github.io/aae.pop/articles/metapopulations.md),
and [multispecies
models](https://aae-stats.github.io/aae.pop/articles/multiple_species.md).

This website includes a [quick start
guide](https://aae-stats.github.io/aae.pop/articles/get_started.md),
examples, vignettes, and [package
documentation](https://aae-stats.github.io/aae.pop/reference/index.md).
Template models have been developed for several species. These are
included in the `aae.pop.templates` package and are described
[here](https://aae-stats.github.io/aae.pop/articles/templates.md).

You can install the current version of the package from Github:

``` r
remotes::install_github("aae-stats/aae.pop")
```

aae.pop is designed as a general package for population models. Several
existing models are included in the aae.pop.templates package, which can
be installed from Github:

``` r
remotes::install_github("aae-stats/aae.pop.templates")
```
