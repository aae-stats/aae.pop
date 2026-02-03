#' @name density_functions
#' @title Common forms of density dependence
#' @description Use pre-defined forms of density dependence
#'   based on common density-dependence functions.
NULL

#' @rdname density_functions
#'
#' @export
#'
#' @param k carrying capacity used to define models of
#'   density dependence. See details for
#'   currently implemented models and their parameters.
#' @param exclude vector of classes to exclude from calculation
#'   of total population density. Defaults to NULL, in which
#'   case all classes are used
#' @param type one of \code{rate} or \code{population}, specifying whether
#'   \code{plot_density_functions} displays the scaling factor for a single
#'   vital rate or the change in population size from one time step to the
#'   next
#' @param \dots ignored
#'
#' @details Additional functions are provided to define common
#'   forms of density dependence. Currently implemented models
#'   are the Ricker model, Beverton-Holt model, a proportional ceiling
#'   model, a theta-logistic model, and alternative parameterisations
#'   of the Ricker and Beverton-Holt models. The functions listed here take
#'   an input \code{k} (the carrying capacity) and an optional parameter
#'   \code{exclude} specifying any classes to exclude from the summed
#'   population size used in density-dependence calculations. However,
#'   these functions return functions, which take a different set of
#'   parameters (see \code{returns}, below).
#'
#'   The \code{plot_density_functions} function produces a basic plot
#'   illustrating the effects of each density dependence model on vital
#'   rates.
#'
#' @returns functions that can be used with \code{\link{density_dependence}}
#'   to specify common models of density dependence. The carrying capacity
#'   \code{k} is specified in the initial function call, but density functions
#'   also take parameters specific to each model.
#'
#'   The \code{ceiling} model function takes no additional parameters.
#'
#'   The \code{ricker} and \code{beverton_holt} model functions each take a
#'   parameter \code{theta}, used to scale the population size in density
#'   calculations (or, alternatively, to inversely scale the carrying
#'   capacity). This parameter is largely redundant but can be used to adjust
#'   \code{k} as an argument to \code{\link{simulate}} rather than define a
#'   new density dependence model.
#'
#'   The \code{ricker2} model function takes a parameter \code{r},
#'   specifying the magnitude of effect when the population is at zero.
#'   This parameter defaults to a value of 1, which multiplies vital
#'   rates by \code{exp(1)} when the population is at zero.
#'
#'   The \code{beverton_holt2} model function takes a parameter \code{alpha},
#'   specifying the magnitude of effect when the population is at carrying
#'   capacity. This parameter defaults to a value of 2, which sets vital
#'   rates to their baseline value when at carrying capacity.
#'
#'   The \code{theta_logistic} model function takes up to three parameters:
#'   \code{s0}, \code{sk}, and (optionally) \code{theta}. The parameter
#'   \code{s0} specifies the scaling of the vital rate when the population
#'   is at zero, \code{sk} specifies this scaling at carrying capacity, and
#'   \code{theta} specifies the strength of density dependence. The default
#'   value of \code{theta} is 1, which is equivalent to a Beverton-Holt
#'   model (with intercepts determined by \code{s0} and \code{sk}). Increasing
#'   \code{theta} generates larger declines in vital rates as a population
#'   increases.
#'
#' @examples
#' # define a population matrix (columns move to rows)
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- dynamics(popmat)
#'
#' # add some density dependence
#' dd <- density_dependence(
#'   masks = reproduction(popmat, dims = 4:5),
#'   funs = ricker(1000)
#' )
#'
#' # update the dynamics object
#' dyn <- update(dyn, dd)
#'
#' # simulate trajectories
#' sims <- simulate(dyn, nsim = 100, options = list(ntime = 50))
#'
#' # and plot
#' plot(sims)
#'
#' # replace with a theta-logistic model
#' dd <- density_dependence(
#'   masks = reproduction(popmat, dims = 4:5),
#'   funs = theta_logistic(1000)
#' )
#'
#' # update the dynamics object
#' dyn <- update(dyn, dd)
#'
#' # simulate trajectories with args passed to theta_logistic
#' sims <- simulate(
#'   dyn,
#'   nsim = 100,
#'   args = list(density_dependence = list(s0 = 1.5, sk = 1, theta = 1.5)),
#'   options = list(ntime = 50)
#' )
#'
#' # and plot
#' plot(sims)
beverton_holt <- function(k, exclude = NULL) {
  function(x, n, theta = 1, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    x / (1 + theta * x * sum(n) / k)
  }
}

#' @rdname density_functions
#'
#' @export
ricker <- function(k, exclude = NULL) {
  function(x, n, theta = 1, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    x * exp(1 - theta * sum(n) / k) / exp(1)
  }
}

#' @rdname density_functions
#'
#' @export
theta_logistic <- function(k, exclude = NULL) {
  function(x, n, s0, sk, theta = 1, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    a <- (s0 / sk) - 1
    x * (s0 / (1 + a * (sum(n) / k) ^ theta))
  }
}

#' @rdname density_functions
#'
#' @export
ceiling_density <- function(k, exclude = NULL) {
  function(x, n, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    x * pmin(k / sum(n), 1)
  }
}

#' @rdname density_functions
#'
#' @export
beverton_holt2 <- function(k, exclude = NULL) {
  function(x, n, alpha = 2, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    x * (alpha / (1 + sum(n) / k))
  }
}

#' @rdname density_functions
#'
#' @export
ricker2 <- function(k, exclude = NULL) {
  function(x, n, r = 1, ...) {
    if (!is.null(exclude)) {
      n <- n[-exclude]
    }
    x * exp(r * (1 - sum(n) / k))
  }
}

#' @rdname density_functions
#'
#' @export
plot_density_functions <- function(type = "rate", ...) {

  # check type
  if (!type %in% c("rate", "population")) {
    stop("type must be one of `rate` or `population`", call. = FALSE)
  }

  # set up a sequence of abundances
  x <- seq(0, 200, length = 100)
  k <- 100

  # calculate the density effect under each model
  vals <- list(
    density_vector(x, 1, beverton_holt(k = k)),
    density_vector(x, 1, beverton_holt2(k = k)),
    density_vector(x, 1, ricker(k = k)),
    density_vector(x, 1, ricker2(k = k)),
    density_vector(x, 1, theta_logistic(k = k), s0 = 1.8, sk = 1, theta = 1),
    density_vector(x, 1, theta_logistic(k = k), s0 = 1.8, sk = 1, theta = 3),
    density_vector(x, 1, ceiling_density(k = k)),
    rep(1, 100)
  )

  # rescale if type is set to population
  xlab <- "Population abundance"
  ylab <- "Effect on vital rate"
  ylim = c(0, exp(1))
  legend_pos <- "topright"
  if (type == "population") {
    vals <- lapply(vals, \(.x, .y) .y * .x, .y = x)
    xlab <- "Population abundance (time t)"
    ylab <- "Population abundance (time t + 1)"
    ylim <- c(0, max(x))
    legend_pos <- "topleft"
  }

  # plot these as basic line plots
  col_pal <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  plot(
    vals[[1]] ~ x,
    type = "n",
    col = col_pal[1],
    ylim = ylim,
    las = 1,
    bty = "l",
    xlab = xlab,
    ylab = ylab
  )

  # add each curve one-by-one
  for (i in seq_along(vals))
    lines(vals[[i]] ~ x, col = col_pal[i], lwd = 2.5)

  # add a legend
  legend(
    legend_pos,
    legend = c(
      "beverton_holt",
      "beverton_holt2",
      "ricker",
      "ricker2",
      "theta_logistic (theta = 1)",
      "theta_logistic (theta = 3)",
      "ceiling_density",
      "none"
    ),
    col = col_pal,
    lwd = 2.5,
    lty = 1,
    bty = "n",
    inset = c(0, 0)
  )

  # return nothing (silently)
  invisible(NULL)

}

# internal function: vectorise density calculations
density_vector <- function(n, x, fn, ...) {
  sapply(n, \(.x) fn(x, .x, ...))
}

