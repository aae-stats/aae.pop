## aae.pop: population dynamics models in R

aae.pop is a plug-and-play tool to simulate matrix population models. It's designed to be fast, flexible, and easily adapted to different model structures, including [complex demographic processes](articles/including_processes.html), [metapopulations](articles/metapopulations.html), and [multispecies models](articles/multiple_species.html).

This website includes a [quick start guide](articles/get_started.html), examples, vignettes, and [package documentation](reference/index.html). Template models have been developed for several species. These are included in the `aae.pop.templates` package and are described [here](articles/templates.html).

You can install the current version of the package from CRAN 

``` r 
install.packages("aae.pop")
```

The latest, development version of `aae.pop` can be installed from Github with the `remotes` package:

``` r
remotes::install_github("aae-stats/aae.pop")
```

aae.pop is designed as a general package for population models. Several existing models are included in the aae.pop.templates package, which can be installed from Github:

``` r
remotes::install_github("aae-stats/aae.pop.templates")
```


<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/aae.pop)](https://CRAN.R-project.org/package=aae.pop)
![R-CMD-check](https://github.com/aae-stats/aae.pop/actions/workflows/check-standard.yaml/badge.svg)
[![Codecov test coverage](https://codecov.io/github/aae-stats/aae.pop/main/graph/badge.svg)](https://app.codecov.io/github/aae-stats/aae.pop/)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
<!-- badges: end -->
