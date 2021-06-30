# functions to summarise simulated population objects

#' @name pr_extinct
#'
#' @title Calculate (quasi-)extinction risk for a \code{\link{simulate}}
#'   object
#'
#' @export
#'
#' @param sims an object returned from \code{\link{simulate}}
#' @param threshold \code{integer} or \code{numeric} denoting the
#'   threshold population size below which a population is considered
#'   functionally extinct. Defaults to \code{0}
#' @param subset \code{integer} vector denoting the population classes
#'   to include in calculation of population abundance. Defaults to
#'   all classes
#' @param times \code{integer} vector specifying generations to
#'   include in calculation of extinction risk. Defaults to all
#'   simulated generations
#'
#' @details Quasi-extinction risk is the probability of decline
#'   below some specified abundance threshold. This probability is
#'   extracted from a \code{\link{simulate}} object as the
#'   proportion of replicate trajectories that fall below this
#'   threshold at any time step within a set period. Abundances
#'   can be specified for all population classes or for a subset
#'   of classes.
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
#' sims <- simulate(dyn, nsim = 1000)
#'
#' # calculate quasi-extinction risk at a threshold population size
#' #   of 100 individuals
#' pr_extinct(sims, threshold = 100)
#'
#' # repeat previous but focused on 4 and 5 year olds only
#' pr_extinct(sims, threshold = 100, subset = 4:5)
#'
#' # repeat previous but ignore first 10 years
#' pr_extinct(sims, threshold = 100, times = 11:51)
pr_extinct <- function(sims, threshold = 0, subset = NULL, times = NULL) {

  # check input object
  if (!"simulation" %in% class(sims))
    stop("pr_extinct is defined for simulation objects", call. = FALSE)

  # check subset
  if (!is.null(subset))
    sims <- subset(sims, subset = subset)

  # check times
  if (!is.null(times))
    sims <- sims[, , times, drop = FALSE]

  # calculate  total abundance
  abund <- apply(sims, c(1, 3), sum)

  # calculate proportion of trajectories that fall below threshold
  out <- mean(
    apply(abund, 1, function(x, y) any(x < y), y = threshold)
  )

  # return
  out

}

#' @name risk_curve
#'
#' @title Calculate (quasi-)extinction risk at multiple thresholds
#'   for a \code{\link{simulate}} object
#'
#' @export
#'
#' @param sims an object returned from \code{\link{simulate}}
#' @param threshold \code{integer} or \code{numeric} vector
#'   denoting the set of threshold population sizes used to define
#'   the risk curve. Defaults to \code{n} evenly spaced values
#'   from 0 to the maximum observed abundance
#' @param subset \code{integer} vector denoting the population classes
#'   to include in calculation of population abundance. Defaults to
#'   all classes
#' @param times \code{integer} vector specifying generations to
#'   include in calculation of extinction risk. Defaults to all
#'   simulated generations
#' @param n \code{integer} specifying number of threshold values
#'   to use in default case when \code{threshold} is not specified.
#'   Defaults to 100
#'
#' @details Risk curves represent \code{\link{pr_extinct}} at multiple
#'   threshold population sizes simultaneously. This gives an expression
#'   of risk of population declines below a range of values. Risk curves
#'   are extracted from a \code{\link{simulate}} object as the
#'   proportion of replicate trajectories that fall below each
#'   threshold value at any time step within a set period. Abundances
#'   can be specified for all population classes or for a subset
#'   of classes.
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
#' sims <- simulate(dyn, nsim = 1000)
#'
#' # calculate risk curve
#' risk_curve(sims)
#'
#' # calculate risk curve for multiple thresholds between 0 and 100
#' risk_curve(sims, threshold = seq(0, 100, by = 1))
#'
#' # calculate risk curve for 4 and 5 year olds only
#' risk_curve(sims, subset = 4:5)
#'
#' # calculate risk curve but ignore first 10 years
#' risk_curve(sims, times = 11:51)
risk_curve <- function(
  sims, threshold = NULL, subset = NULL, times = NULL, n = 100
) {

  # check input object
  if (!"simulation" %in% class(sims))
    stop("risk_curve is defined for simulation objects", call. = FALSE)

  # check subset
  if (!is.null(subset))
    sims <- subset(sims, subset = subset)

  # check times
  if (!is.null(times))
    sims <- sims[, , times, drop = FALSE]

  # set default threshold if required
  if (is.null(threshold)) {

    # calculate  total abundance
    abund <- apply(sims, c(1, 3), sum)

    # calculate maximum abundance of any replicate
    max_n <- max(apply(abund, 1, max))

    # define a threshold as n evenly spaced values
    #   between 0 and the max abundance
    threshold <- seq(0, max_n, length = n)

  }

  # calculate risk at each threshold, noting that
  #   subset and times are accounted for above
  out <- sapply(
    threshold,
    function(x, y) pr_extinct(sims = y, threshold = x),
    y = sims
  )

  # add names to output so thresholds are known
  names(out) <- round(threshold, digits = 0)

  # return
  out

}

#' @name emps
#'
#' @title Calculate expected minimum population size (EMPS)
#'   for a \code{\link{simulate}} object
#'
#' @export
#'
#' @param sims an object returned from \code{\link{simulate}}
#' @param subset \code{integer} vector denoting the population classes
#'   to include in calculation of population abundance. Defaults to
#'   all classes
#' @param times \code{integer} vector specifying generations to
#'   include in calculation of extinction risk. Defaults to all
#'   simulated generations
#' @param fun \code{function} used to calculate average over all
#'   replicate trajectories. Defaults to \code{mean}. Alternatives
#'   might include \code{median} or \code{min}
#' @param \dots additional arguments passed to \code{fun}
#'
#' @details Expected minimum population size (EMPS) is the average
#'   minimum value of all replicate trajectories. This value represents
#'   an expected lower bound on population sizes over all generations,
#'   accounting for variation among replicates. Abundances
#'   can be specified for all population classes or for a subset
#'   of classes.
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
#' sims <- simulate(dyn, nsim = 1000)
#'
#' # calculate expected minimum population size
#' emps(sims)
#'
#' # calculate expected minimum population size for 4 and 5 year
#' #   olds only
#' emps(sims, subset = 4:5)
#'
#' # calculate expected minimum population size but ignore first 10 years
#' emps(sims, times = 11:51)
#'
#' # calculate expected minimum population size based on median
#' emps(sims, fun = median)
emps <- function(sims, subset = NULL, times = NULL, fun = mean, ...) {

  # check input object
  if (!"simulation" %in% class(sims))
    stop("emps is defined for simulation objects", call. = FALSE)

  # check subset
  if (!is.null(subset))
    sims <- subset(sims, subset = subset)

  # check times
  if (!is.null(times))
    sims <- sims[, , times, drop = FALSE]

  # calculate  total abundance
  abund <- apply(sims, c(1, 3), sum)

  # calculate minimum abundance in each trajectory
  min_n <- apply(abund, 1, min)

  # calculate and return average minimum over all trajectories
  fun(min_n, ...)

}

#' @name exps
#'
#' @title Calculate expected population size for a
#'   \code{\link{simulate}} object based on
#'   generic functions (ExPS)
#'
#' @export
#'
#' @param sims an object returned from \code{\link{simulate}}
#' @param subset \code{integer} vector denoting the population classes
#'   to include in calculation of population abundance. Defaults to
#'   all classes
#' @param times \code{integer} vector specifying generations to
#'   include in calculation of extinction risk. Defaults to all
#'   simulated generations
#' @param fun_within \code{function} used to summarise a single
#'   trajectory. Must return a single value. Defaults to \code{mean}
#' @param fun_among \code{function} used to summarise over all
#'   replicate trajectories. Defaults to \code{mean}. Alternatives
#'   might include \code{median} or \code{min}
#' @param \dots additional arguments passed to \code{fun_within} and
#'   \code{fun_among}. If these conflict, a wrapper function could
#'   be used to define expected arguments for each function
#'
#' @details Expected population size (ExPS) is a highly
#'   flexible generalisation of \code{\link{emps}} and
#'   represents a two-level summary that first summarises
#'   individual population trajectories and then summarisies
#'   these values over all replicates. Abundances
#'   can be specified for all population classes or for a subset
#'   of classes.
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
#' sims <- simulate(dyn, nsim = 1000)
#'
#' # calculate expected population size
#' exps(sims)
#'
#' # calculate expected population size for 4 and 5 year
#' #   olds only
#' exps(sims, subset = 4:5)
#'
#' # calculate expected population size but ignore first 10 years
#' exps(sims, times = 11:51)
#'
#' # calculate expected population size based on median
#' exps(sims, fun_among = median)
#'
#' # calculate expected maximum population size based on median
#' exps(sims, fun_within = max, fun_among = median)
#'
#' # calculate exps with conflicting quantile functions, handling
#' #   conflicting arguments with wrapper functions
#' quant1 <- function(x, p1, ...) {quantile(x, prob = p1)}
#' quant2 <- function(x, p2, ...) {quantile(x, prob = p2)}
#' exps(
#'   sims, fun_within = quant1, fun_among = quant2, p1 = 0.25, p2 = 0.75
#' )
exps <- function(
  sims, subset = NULL, times = NULL, fun_within = mean, fun_among = mean, ...
) {

  # check input object
  if (!"simulation" %in% class(sims))
    stop("exps is defined for simulation objects", call. = FALSE)

  # check subset
  if (!is.null(subset))
    sims <- subset(sims, subset = subset)

  # check times
  if (!is.null(times))
    sims <- sims[, , times, drop = FALSE]

  # calculate  total abundance
  abund <- apply(sims, c(1, 3), sum)

  # calculate summary function within in each trajectory
  summary_n <- apply(abund, 1, fun_within, ...)

  # calculate and return average minimum over all trajectories
  fun_among(summary_n, ...)

}
