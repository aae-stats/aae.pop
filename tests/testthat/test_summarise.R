context("summarise")

# setup: simulate some data to test with
nstage <- 5
mat <- matrix(0, nrow = nstage, ncol = nstage)
mat[reproduction(mat, dims = 4:5)] <- rpois(2, 20)
mat[survival(mat)] <- plogis(rnorm(nstage))
mat[transition(mat)] <- plogis(rnorm(nstage - 1))

# add covariate effects
ntime <- 35
xsim <- rnorm(ntime)
cov_fn <- function(mat, x) {
  mat * plogis(x)
}
cov_eff <- covariates(
  masks = survival(mat),
  funs = cov_fn
)
rep_cov <- replicated_covariates(
  masks = survival(mat),
  funs = cov_fn
)

# add density dependence
dd_masks <- list(reproduction(mat, dims = 4:5))
dd_fns <- list(
  function(x, n) x * (1 / (1 + 1e-4 * sum(n)))
)
dd <- density_dependence(masks = dd_masks, funs = dd_fns)

# simulate from basic model
nsim <- 100
dyn <- dynamics(mat, cov_eff, dd)
init_set <- matrix(rpois(nstage * nsim, lambda = 5), ncol = nstage)
sims <- simulate(
  dyn,
  nsim = nsim,
  init = init_set,
  options = list(ntime = ntime, tidy_abundances = floor),
  args = list(
    covariates = format_covariates(xsim)
  )
)

# calculate summed abundances for use in targets
abund <- apply(sims, c(1, 3), sum)
abund_sub <- apply(sims[, 3:5, ], c(1, 3), sum)

test_that("pr_extinct works with multiple settings", {

  # defaults
  value <- pr_extinct(sims = sims)
  target <- mean(apply(abund, 1, function(x) any(x < 0)))
  expect_equal(value, target)

  # increase threshold
  value <- pr_extinct(sims = sims, threshold = 100)
  target <- mean(apply(abund, 1, function(x) any(x < 100)))
  expect_equal(value, target)

  # negative threshold
  value <- pr_extinct(sims = sims, threshold = -100)
  target <- mean(apply(abund, 1, function(x) any(x < -100)))
  expect_equal(value, target)

  # subset of classes
  value <- pr_extinct(sims = sims, subset = 3:5)
  target <- mean(apply(abund_sub, 1, function(x) any(x < 0)))
  expect_equal(value, target)

  # subset of times
  value <- pr_extinct(sims = sims, times = 20:35)
  target <- mean(apply(abund[, 20:35], 1, function(x) any(x < 0)))
  expect_equal(value, target)

  # combine all
  value <- pr_extinct(sims = sims, threshold = 42000, subset = 3:5, times = 20:35)
  target <- mean(apply(abund_sub[, 20:35], 1, function(x) any(x < 42000)))
  expect_equal(value, target)

})

test_that("risk_curve works with multiple settings", {

  # defaults
  value <- risk_curve(sims = sims)
  risk_seq <- seq(0, max(abund), length = 100)
  target <- sapply(
    risk_seq, function(x, y) pr_extinct(y, threshold = x), y = sims
  )
  names(target) <- round(risk_seq, 0)
  expect_equal(value, target)

  # change threshold
  value <- risk_curve(sims = sims, threshold = c(5, 10, 20))
  risk_seq <- c(5, 10, 20)
  target <- sapply(
    risk_seq, function(x, y) pr_extinct(y, threshold = x), y = sims
  )
  names(target) <- round(risk_seq, 0)
  expect_equal(value, target)

  # subset of classes
  value <- risk_curve(sims = sims, subset = 3:5)
  risk_seq <- seq(0, max(abund_sub), length = 100)
  target <- sapply(
    risk_seq,
    function(x, y) pr_extinct(y, subset = 3:5, threshold = x),
    y = sims
  )
  names(target) <- round(risk_seq, 0)
  expect_equal(value, target)

  # subset of times
  value <- risk_curve(sims = sims, times = 20:35)
  risk_seq <- seq(0, max(abund[, 20:35]), length = 100)
  target <- sapply(
    risk_seq,
    function(x, y) pr_extinct(y, times = 20:35, threshold = x),
    y = sims
  )
  names(target) <- round(risk_seq, 0)
  expect_equal(value, target)

  # combine all
  value <- risk_curve(sims = sims, subset = 3:5, times = 20:35, n = 45)
  risk_seq <- seq(0, max(abund_sub[, 20:35]), length = 45)
  target <- sapply(
    risk_seq,
    function(x, y) {
      pr_extinct(y, subset = 3:5, times = 20:35, threshold = x)
    },
    y = sims
  )
  names(target) <- round(risk_seq, 0)
  expect_equal(value, target)

})

test_that("get_pdf works with multiple settings", {

  # defaults
  value <- get_pdf(sims)
  target <- density(apply(abund, 1, min), n = 100)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # subset of classes
  value <- get_pdf(sims, subset = 3:5)
  target <- density(apply(abund_sub, 1, min), n = 100)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # subset of times
  value <- get_pdf(sims, times = 20:35)
  target <- density(apply(abund[, 20:35], 1, min), n = 100)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # change number of bins
  value <- get_pdf(sims, n = 20)
  target <- density(apply(abund, 1, min), n = 20)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # change function
  value <- get_pdf(sims, fn = max)
  target <- density(apply(abund, 1, max), n = 100)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # args to function
  value <- get_pdf(sims, fn = quantile, probs = 0.25)
  target <- density(apply(abund, 1, quantile, probs = 0.25), n = 100)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

  # combine all
  value <- get_pdf(sims, subset = 3:5, times = 20:35, fn = max, n = 30)
  target <- density(apply(abund_sub[, 20:35], 1, max), n = 30)
  target <- data.frame(prob = target$y, value = target$x)
  expect_equal(value, target)

})

test_that("get_cdf works with multiple settings", {

  # defaults
  value <- get_cdf(sims)
  cdf_seq <- seq(0, 1, length = 100 + 1)
  target <- apply(abund, 1, min)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # subset of classes
  value <- get_cdf(sims, subset = 3:5)
  cdf_seq <- seq(0, 1, length = 100 + 1)
  target <- apply(abund_sub, 1, min)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # subset of times
  value <- get_cdf(sims, times = 20:35)
  cdf_seq <- seq(0, 1, length = 100 + 1)
  target <- apply(abund[, 20:35], 1, min)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # change number of bins
  value <- get_cdf(sims, n = 20)
  cdf_seq <- seq(0, 1, length = 20 + 1)
  target <- apply(abund, 1, min)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # change function
  value <- get_cdf(sims, fn = max)
  cdf_seq <- seq(0, 1, length = 100 + 1)
  target <- apply(abund, 1, max)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # arguments to function
  value <- get_cdf(sims, fn = quantile, probs = 0.25)
  cdf_seq <- seq(0, 1, length = 100 + 1)
  target <- apply(abund, 1, quantile, probs = 0.25)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

  # combine all
  value <- get_cdf(sims, subset = 3:5, times = 20:35, fn = max, n = 30)
  cdf_seq <- seq(0, 1, length = 30 + 1)
  target <- apply(abund_sub[, 20:35], 1, max)
  target <- quantile(target, probs = cdf_seq)
  names(target) <- NULL
  target <- data.frame(prob = cdf_seq, value = target)
  expect_equal(value, target)

})

test_that("emps works with multiple settings", {

  # defaults
  value <- emps(sims)
  target <- mean(apply(abund, 1, min))
  expect_equal(value, target)

  # subset of classes
  value <- emps(sims, subset = 3:5)
  target <- mean(apply(abund_sub, 1, min))
  expect_equal(value, target)

  # subset of times
  value <- emps(sims, times = 20:35)
  target <- mean(apply(abund[, 20:35], 1, min))
  expect_equal(value, target)

  # change function
  value <- emps(sims, fun = max)
  target <- max(apply(abund, 1, min))
  expect_equal(value, target)

  # arguments to function
  value <- emps(sims, fun = quantile, probs = 0.25)
  target <- quantile(apply(abund, 1, min), probs = 0.25)
  expect_equal(value, target)

  # combine all
  value <- emps(sims, subset = 3:5, times = 20:35, fun = max)
  target <- max(apply(abund_sub[, 20:35], 1, min))
  expect_equal(value, target)

})

test_that("exps works with multiple settings", {

  # defaults
  value <- exps(sims)
  target <- mean(apply(abund, 1, mean))
  expect_equal(value, target)

  # subset of classes
  value <- exps(sims, subset = 3:5)
  target <- mean(apply(abund_sub, 1, mean))
  expect_equal(value, target)

  # subset of times
  value <- exps(sims, times = 20:35)
  target <- mean(apply(abund[, 20:35], 1, mean))
  expect_equal(value, target)

  # change among function
  value <- exps(sims, fun_among = max)
  target <- max(apply(abund, 1, mean))
  expect_equal(value, target)

  # change within function
  value <- exps(sims, fun_within = max)
  target <- mean(apply(abund, 1, max))
  expect_equal(value, target)

  # arguments to within function
  value <- exps(sims, fun_within = quantile, probs = 0.25, fun_among = max)
  target <- max(apply(abund, 1, quantile, probs = 0.25))
  expect_equal(value, target)

  # arguments to among function
  value <- exps(sims, fun_among = quantile, probs = 0.25, fun_within = max)
  target <- quantile(apply(abund, 1, max), probs = 0.25)
  expect_equal(value, target)

  # combine all
  value <- exps(
    sims,
    subset = 3:5,
    times = 20:35,
    fun_among = quantile,
    probs = 0.25,
    fun_within = max
  )
  target <- quantile(apply(abund_sub[, 20:35], 1, max), probs = 0.25)
  expect_equal(value, target)

  # arguments to both functions
  fn1 <- function(x, p1, ...) quantile(x, probs = p1)
  fn2 <- function(x, p2, ...) quantile(x, probs = p2)
  value <- exps(sims, fun_among = fn1, p1 = 0.25, fun_within = fn2, p2 = 0.5)
  target <- quantile(apply(abund, 1, quantile, probs = 0.5), probs = 0.25)
  expect_equal(value, target)

})

test_that("all summarise methods error informatively", {

  # pr_extinct
  expect_error(
    pr_extinct(rnorm(10)),
    "is defined for simulation objects"
  )

  # risk_curve
  expect_error(
    risk_curve(rnorm(10)),
    "is defined for simulation objects"
  )

  # get_pdf
  expect_error(
    get_pdf(rnorm(10)),
    "is defined for simulation objects"
  )

  # get_cdf
  expect_error(
    get_cdf(rnorm(10)),
    "is defined for simulation objects"
  )

  # emps
  expect_error(
    emps(rnorm(10)),
    "is defined for simulation objects"
  )

  # exps
  expect_error(
    exps(rnorm(10)),
    "is defined for simulation objects"
  )

})
