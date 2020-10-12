#' @name updaters
#' @title Functions for a single time-step update (matrix %*% vector)
#' @description Define how population abundances are updated from
#'   one time step to the next. Functions can take any form but will
#'   only be vectorised across replicates in limited situations.
NULL

#' @rdname templates
#'
#' @export
#'
#' @param pop current state of the population
#' @param mat matrix of vital rates used to update population state
#'
#' @details Update abundances with a direct matrix multiplication
#'   that does not include any form of demographic stochasticity.
#'   This is the fastest update option and will vectorise across
#'   replicates if the population matrix is not expanded by
#'   \code{\link{environmental_stochasticity}} or
#'   \code{\link{density_dependence}}.
update_crossprod <- function(pop, mat) {
  tcrossprod(pop, mat)
}

#' @rdname templates
#'
#' @export
#'
#' @details Update abundances with a direct RNG draw
#   that combines update with demographic stochasticity,
#   assuming a Leslie matrix.
update_binomial_leslie <- function(pop, mat) {

  if (!all(pop%%1 == 0)) {
    stop("some abundances are not integers, so cannot be used ",
         "with update_binomial_leslie. Check options()$tidy_abundances ",
         "and update with an appropriate method (e.g. floor)",
         call. = FALSE)
  }

  pop_nm1 <- pop[-length(pop)]

  vals <- tcrossprod(pop, mat)
  probs <- vals[-1] / pop_nm1
  probs[pop_nm1 == 0] <- 0

  c(
    rpois(1, lambda = vals[1]),
    rbinom(length(probs), size = pop_nm1, prob = probs)
  )

}

#' @rdname templates
#'
#' @importFrom mc2d rmultinomial
#'
#' @export
#'
#' @details Update abundances with a direct RNG draw
#'   that combines update with demographic stochasticity,
#'   allowing for general matrix forms (slower than
#'   \code{update_binomial_leslie}).
update_multinomial <- function(pop, mat) {

  if (!all(pop%%1 == 0)) {
    stop("some abundances are not integers, so cannot be used ",
         "with update_multinomial. Check options()$tidy_abundances ",
         "and update with an appropriate method (e.g. floor)",
         call. = FALSE)
  }

  n <- length(pop)

  recruits <- mat[1, ]
  recruits[1] <- 0
  recruits <- rpois(1, lambda = crossprod(recruits, pop))

  mat[reproduction(mat)] <- 0

  # add class for dead individuals
  mat <- rbind(mat, 1 - colSums(mat))

  out <- matrix(0, nrow = n, ncol = n + 1)

  out <- mc2d::rmultinomial(n, size = pop, prob = t(mat))
  out <- apply(out[, 1:n], 2, sum)

  out[1] <- out[1] + recruits

  out

}
