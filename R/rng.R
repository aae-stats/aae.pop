#' @name rng
#' @title Random number generators not available in existing R
#'   packages
#' @description Draw random numbers from unusual distributions,
#'   such as on the unit or non-negative real line with known
#'   means and standard deviations.
NULL

#' @rdname rng
#'
#' @importFrom nleqslv nleqslv
#' @importFrom cubature hcubature
#' @importFrom stats pnorm rnorm qnorm dnorm cov2cor uniroot integrate
#'
#' @export
#'
#' @param n number of random draws to simulate. Each draw is a vector of
#'   values with length equal to \code{length(mean)} and
#'   \code{length(sd)} and the resulting output has \code{n} rows
#'   and \code{length(mean)} columns
#' @param mean vector of mean values on the unit scale
#' @param sd vector of positive standard deviations
#' @param Sigma optional covariance matrix with dimensions of
#'   \code{length(mean)} by \code{length(mean)} defining covariances
#'   between each pair of values in \code{mean}. Note that only
#'   the correlation structure is retained from \code{Sigma},
#'   so that standard deviations are still required
#' @param Omega optional correlation matrix with dimensions of
#'   \code{length(mean)} by \code{length(mean)} defining correlations
#'   between each pair of values in \code{mean}
#' @param perfect_correlation \code{logical}, if \code{TRUE}
#'   and \code{Sigma} and \code{Omega} are \code{NULL}, then all
#'   values in each replicate (row) are perfectly correlated with known
#'   mean and standard deviation. If \code{FALSE}, then all values
#'   in each replicate are completely uncorrelated
#'
#' @details Function to simulate values on unit interval with fixed
#'   mean, sd, and correlation
rmultiunit <- function(n, mean, sd, Sigma = NULL, Omega = NULL, perfect_correlation = FALSE) {

  # how many parameters are we dealing with?
  npar <- length(mean)

  # convert covariance to correlation if provided
  if (!is.null(Sigma) & is.null(Omega))
    Omega <- cov2cor(Sigma)

  # calculate mean and sd on the real line
  real_params <- unit_to_real(unit_mean = mean,
                              unit_sd = sd)

  # calculate covariance if correlation/covariance matrix provided
  if (!is.null(Omega)) {

    # estimate full covariance matrix on real line
    Sigma <- unit_to_real_covar(corr = Omega,
                                unit_mean = mean,
                                unit_sd = sd,
                                real_params = real_params)

    # take Cholesky decomposition to get a triangular matrix
    Sigma_chol <- t(chol(Sigma))

  }

  # simulate random values from a standard normal
  z_variates <- matrix(rnorm(npar * n), ncol = n)

  # simpler return if correlations not required
  if (is.null(Omega)) {

    # do we want perfectly correlated or uncorrelated?
    if (perfect_correlation) {
      out <- t(pnorm(real_params[, 1] + real_params[, 2] %o% z_variates[1, ]))
    } else {
      out <- t(pnorm(real_params[, 1] + sweep(z_variates, 1, real_params[, 2], "*")))
    }

  } else {

    # combine with correlations and means to give full variates
    out <- t(pnorm(real_params[, 1] + Sigma_chol %*% z_variates))

  }

  # return
  out

}

# equation based on the expected mean
phi_mean <- function(x, y) {
  pnorm(x / sqrt(1 + y ^ 2))
}

# equation based on the expected variance
variance_calc <- function(z, x, y) {
  (pnorm(x + y * z) ^ 2) * dnorm(z)
}

# internal function to integrate over a wide range of values
phi_var <- function(x, y, lower = -10, upper = 10) {
  integrate(f = variance_calc, lower = lower, upper = upper, x = x, y = y)$value
}

# define the system of equations to solve
f_xy <- function(x, p, s) {
  c(phi_mean(x[1], x[2]) - p,
    phi_var(x[1], x[2]) - (p * p + s * s))
}

# 2D Newton solver to find the roots of two non-linear equations
solve_nl <- function(x, fn, ...) {

  init <- c(qnorm(x[1]), x[2] / dnorm(qnorm(x[1])))

  nleqslv::nleqslv(x = init,
                   fn = fn,
                   p = x[1],
                   s = x[2],
                   ...)$x

}

# function to calculate mean and sd on real line
unit_to_real <- function(unit_mean, unit_sd) {
  t(apply(cbind(unit_mean, unit_sd), 1, solve_nl, fn = f_xy))
}

# integrand to define correlation coefficients
rho_int <- function(x, mean_i, mean_j, sd_i, sd_j, rho, rho2) {
  matrix(
    pnorm(mean_i + sd_i * x[1, ]) *
      pnorm(mean_j + sd_j * x[2, ]) *
      (1 / (2 * pi * sqrt(1 - rho2))) *
      exp(-(1 / (2 * (1 - rho2))) * (x[1, ] ^ 2 - 2 * rho * x[1, ] * x[2, ] + x[2, ] ^ 2)),
    ncol = ncol(x)
  )
}

# equation to solve to find correlation coefficient
f_x <- function(x, fn, arg) {

  # integrate over wide range of values in both dimensions
  int_value <- cubature::hcubature(
    fn,
    lowerLimit = rep(-10, 2),
    upperLimit = rep(10, 2),
    mean_i = arg["mean_i"],
    mean_j = arg["mean_j"],
    sd_i = arg["sd_i"],
    sd_j = arg["sd_j"],
    rho = x,
    rho2 = x ^ 2,
    vectorInterface = TRUE,
    tol = 1e-5
  )

  # return the calculated equation based on the expected correlation
  int_value$integral - (arg["pi"] * arg["pj"] + arg["rij"] * arg["si"] * arg["sj"])

}

# wrap up the root finder into a tidier function
root_finder <- function(arg, f, fn, interval = c(-0.05, 0.05), extend = "yes") {
  uniroot(f = f, interval = interval, extendInt = extend, fn = fn, arg = arg)$root
}

# function to calculate covariance on real line
unit_to_real_covar <- function(corr, unit_mean, unit_sd, real_params) {

  # unpack the inputs so we can vectorise this
  idz <- upper.tri(corr)
  idx <- row(corr)[idz]
  idy <- col(corr)[idz]

  # collate it all into a matrix
  args <- cbind(
    pi = unit_mean[idx],
    pj = unit_mean[idy],
    si = unit_sd[idx],
    sj = unit_sd[idy],
    rij = corr[upper.tri(corr)],
    mean_i = real_params[idx, 1],
    mean_j = real_params[idy, 1],
    sd_i = real_params[idx, 2],
    sd_j = real_params[idy, 2]
  )

  # calculate pairwise correlations
  corr[upper.tri(corr)] <- apply(args, 1, root_finder, f = f_x, fn = rho_int)

  # calculate and return covariance matrix
  corr[lower.tri(corr)] <- corr[upper.tri(corr)]
  diag(corr) <- 1
  diag(real_params[, 2]) %*% corr %*% diag(real_params[, 2])

}
