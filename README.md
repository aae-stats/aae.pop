## aae.pop: simulating multispecies population dynamics in R

## Installing the package

You can install the `aae.pop` package from GitHub. To install from GitHub, you'll need to install the `remotes` R package and use the following lines of code:

```
# install the remotes package if not already installed
install.packages("remotes")

# install the aae.pop package from GitHub
remotes::install_github("aae.stats/aae.pop@master")
```

Once completed, you should be able to load the `aae.pop` package with `library(aae.pop)`.

## Building a basic population model

```
# load the aae.pop package
library(aae.pop)

# create a basic Leslie matrix for a population with
#   five age classes
# the aae.pop package assumes "columns move to rows"
popmat <- rbind(
  c(0,    0,    2,    4,    7),  # reproduction from 3-5 year olds
  c(0.25, 0,    0,    0,    0),  # survival from age 1 to 2
  c(0,    0.45, 0,    0,    0),  # survival from age 2 to 3
  c(0,    0,    0.70, 0,    0),  # survival from age 3 to 4
  c(0,    0,    0,    0.85, 0)   # survival from age 4 to 5
)

# create a population dynamics object with this matrix
#   without any additional processes
popdyn <- dynamics(popmat)

# can plot the model structure if the DiagrammeR package is installed
plot(popdyn)

# simulate from this model with default (random) initial conditions
sims <- simulate(popdyn)

# can plot this
plot(sims, xlab = "Generation", ylab = "Abundance")
```

The vignettes contain several more detailed examples.

Please leave feedback, bug reports or feature requests at the GitHub [issues page](https://github.com/aae-stats/aae.pop/issues). 

[![build status](https://travis-ci.org/aae-stats/aae.pop.svg?branch=master)](https://travis-ci.org/aae-stats/aae.pop) [![codecov.io](https://codecov.io/github/aae-stats/aae.pop/coverage.svg?branch=master)](https://codecov.io/github/aae-stats/aae.pop?branch=master) [![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
