% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density_dependence.R
\name{density_dependence}
\alias{density_dependence}
\alias{density_dependence_n}
\title{Specify density dependence in models of population dynamics}
\usage{
density_dependence(masks, funs, nmask = NULL)

density_dependence_n(masks, funs)
}
\arguments{
\item{masks}{a logical matrix or vector (or list of these)
defining cells affected by \code{funs}. See Details and
\code{\link{masks}}}

\item{funs}{a function or list of functions with one element
for each element of \code{masks}. See Details}

\item{nmask}{logical vector or list of vectors defining
elements of the population vector affected by each
mask-function pair. Intended primarily for internal
use when scaling up processes in
\code{\link{metapopulation}}}
}
\description{
Specify density dependence in vital rates
  (\code{density_dependence}) and in total abundances
  (\code{density_dependence_n}).
}
\details{
\code{density_dependence} specifies standard
  density dependence on vital rates, such as scramble or
  contest competition or allee effects.

  \code{density_dependence_n} is an alternative
  parameterisation of density dependence that acts directly
  on population abundances.

  Masks must be of the same dimension as the population
  dynamics matrix and specify cells influenced by density
  dependence according to \code{funs}. In the case of
  \code{density_dependence_n}, \code{masks} are
  logical vectors with one element for each class.
  Additional details on masks are provided
  in \code{\link{masks}}.

  If using \code{density_depenence}, functions must take at
  least two arguments, a matrix \code{x} and a vector \code{n},
  which represent the population dynamics matrix and the
  population abundances. Functions must return a matrix with
  the same dimensions as \code{x}, modified to reflect the
  effects of current abundances (\code{n}) on
  vital rates.

  In the case of \code{density_dependence_n},
  \code{funs} takes only one argument, the population
  abundances \code{n} following all other updates in a
  given iteration/generation. This allows rescaling of
  population abundances based on total abundance or
  through more complicated functions that depend
  on external arguments (e.g., mass mortality events or
  harvesting).

  Additional arguments to functions are supported and can be
  passed to \code{\link{simulate}} with the \code{args},
  \code{args.dyn}, or \code{args.fn} arguments.
}
\examples{
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
}
