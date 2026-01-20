#' @name updaters
#' @title Functions for a single time-step update
#' @description Define how population abundances are updated from
#'   one time step to the next. Functions can take any form but will
#'   only be vectorised across replicates in limited situations.
NULL

#' @rdname updaters
#'
#' @importFrom mc2d rmultinomial
#'
#' @export
#'
#' @param pop current state of the population
#' @param mat matrix of vital rates used to update population state
#'
#' @details Updaters can be changed through the \code{options}
#'   argument to \code{simulate} and can also be changed
#'   globally for an R session by changing the global option with,
#'   e.g., \code{options(aae.pop_update = update_binomial_leslie)}
#'
#'   \code{update_crossprod} updates abundances with a direct
#'   matrix multiplication that does not include any form of demographic
#'   stochasticity. This is the fastest update option and will vectorise
#'   across replicates if the population matrix is not expanded by
#'   \code{\link{environmental_stochasticity}} or
#'   \code{\link{density_dependence}}.
#'
#' @returns a matrix containing population abundances in each stage of a
#'   matrix population model. Contains one row for each replicate population
#'   trajectory and one column for each population stage
#'
#' @examples
#' # define a basic population
#' nstage <- 5
#' popmat <- matrix(0, nrow = nstage, ncol = nstage)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- dynamics(popmat)
#'
#' # simulate with the default updater
#' sims <- simulate(dyn)
#'
#' # simulate with a multinomial updater
#' sims <- simulate(dyn, options = list(update = update_multinomial))
update_crossprod <- function(pop, mat) {
  tcrossprod(pop, mat)
}

#' @rdname updaters
#'
#' @export
#'
#' @details \code{update_binomial_leslie} updates abundances
#'   with a direct RNG draw that combines update with demographic
#'   stochasticity, assuming a Leslie matrix.
update_binomial_leslie <- function(pop, mat) {
  if (!all((pop %% 1) == 0)) {
    stop("some abundances are not integers, so cannot be used ",
      "with update_binomial_leslie. ",
      "Check options()$aae.pop_tidy_abundances ",
      "and update with an appropriate method (e.g. floor)",
      call. = FALSE
    )
  }

  # need a 1 row matrix if pop is a vector
  if (is.null(dim(pop))) {
    pop <- matrix(pop, nrow = 1)
  }

  # pull out all ages except the last
  pop_nm1 <- pop[, -ncol(pop), drop = FALSE]

  vals <- tcrossprod(pop, mat)
  probs <- vals[, -1] / pop_nm1
  probs[pop_nm1 == 0] <- 0

  cbind(
    rpois(nrow(vals), lambda = vals[, 1]),
    matrix(rbinom(length(probs), size = pop_nm1, prob = probs),
      nrow = nrow(pop_nm1)
    )
  )
}

#' @rdname updaters
#'
#' @importFrom mc2d rmultinomial
#'
#' @export
#'
#' @details \code{update_multinomial} updates abundances with a
#'   direct RNG draw that combines update with demographic stochasticity,
#'   allowing for general matrix forms (slower than
#'   \code{update_binomial_leslie}).
update_multinomial <- function(pop, mat) {
  if (!all((pop %% 1) == 0)) {
    stop("some abundances are not integers, so cannot be used ",
      "with update_multinomial. Check options()$tidy_abundances ",
      "and update with an appropriate method (e.g. floor)",
      call. = FALSE
    )
  }

  # need a 1 row matrix if pop is a vector
  if (is.null(dim(pop))) {
    pop <- matrix(pop, nrow = 1)
  }

  n <- ncol(pop)

  recruits <- mat[1, ]
  recruits[1] <- 0
  recruits <- rpois(nrow(pop), lambda = tcrossprod(recruits, pop))

  mat[reproduction(mat)] <- 0

  # add class for dead individuals
  mat <- rbind(mat, 1 - colSums(mat))

  # correct if total probs exceed 1 and dead probs are negative
  mat[mat < 0] <- 0

  out <- matrix(0, nrow = n, ncol = n + 1)

  out <- t(apply(
    pop,
    1,
    multinomial_internal,
    n = n,
    mat = mat
  ))

  out[, 1] <- out[, 1] + recruits

  out
}

# internal function: single multinomial draw for one replicate
multinomial_internal <- function(n, pop, mat) {
  out <- rmultinomial(n = n, size = pop, prob = t(mat))
  apply(out[, 1:n], 2, sum)
}
