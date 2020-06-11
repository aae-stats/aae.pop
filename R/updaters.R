# functions for a single time-step update (matrix %*% vector)

# internal function: update abundances for one time step
update_crossprod <- function(pop, mat) {
  tcrossprod(pop, mat)
}

# internal function: update abundances with a direct RNG draw
#   that combines update with demographic stochasticity
update_binomial <- function(pop, mat) {

  # can't work out size of rbinom unless pop is an integer
  if (!all(pop %% 1 == 0)) {
    stop("some abundances are not integers, so cannot be used ",
         "with update_binomial. Check options()$tidy_abundances ",
         "and update with an appropriate method (e.g. floor)",
         call. = FALSE)
  }

  # run internal update on each class individually
  pop_tp1 <- apply(mat, 1, update_binomial_internal, n = pop)

  # return outputs
  pop_tp1

}

# internal function: update population abundances with a demographic
#   stochasticity function based on Poisson draws for reproduction and
#   Binomial draws for survival
update_binomial_internal <- function(x, n) {

  # which classes contribute to class i?
  source <- x > 0

  # use Poisson if dealing with creation of new individuals,
  #  Binomial otherwise
  if (all(x < 1)) {
    n_tp1 <- sum(rbinom(sum(source), size = n[source], prob = x[source]))
  } else {
    n_tp1 <- sum(rpois(sum(source), lambda = (n[source] * x[source])))
  }

  # return
  n_tp1

}

