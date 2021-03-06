% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updaters.R
\name{updaters}
\alias{updaters}
\alias{update_crossprod}
\alias{update_binomial_leslie}
\alias{update_multinomial}
\title{Functions for a single time-step update}
\usage{
update_crossprod(pop, mat)

update_binomial_leslie(pop, mat)

update_multinomial(pop, mat)
}
\arguments{
\item{pop}{current state of the population}

\item{mat}{matrix of vital rates used to update population state}
}
\description{
Define how population abundances are updated from
  one time step to the next. Functions can take any form but will
  only be vectorised across replicates in limited situations.
}
\details{
Updaters can be changed through the \code{options}
  argument to \code{simulate} and can also be changed
  globally for an R session by changing the global option with,
  e.g., \code{options(aae.pop_update = update_binomial_leslie)}

  \code{update_crossprod} updates abundances with a direct
  matrix multiplication that does not include any form of demographic
  stochasticity. This is the fastest update option and will vectorise
  across replicates if the population matrix is not expanded by
  \code{\link{environmental_stochasticity}} or
  \code{\link{density_dependence}}.

\code{update_binomial_leslie} updates abundances
  with a direct RNG draw that combines update with demographic
  stochasticity, assuming a Leslie matrix.

\code{update_multinomial} updates abundances with a
  direct RNG draw that combines update with demographic stochasticity,
  allowing for general matrix forms (slower than
  \code{update_binomial_leslie}).
}
\examples{
# define a basic population
nstage <- 5
popmat <- matrix(0, nrow = nstage, ncol = nstage)
popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)

# define a dynamics object
dyn <- dynamics(popmat)

# simulate with the default updater
sims <- simulate(dyn)

# simulate with a multinomial updater
sims <- simulate(dyn, options = list(update = update_multinomial))
}
