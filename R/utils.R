# create a neat list of terms for messages, warnings, and errors
clean_paste <- function(x, final_sep = "and") {

  # want number of terms so we can separate the final term
  nterm <- length(x)

  # only clean if x has length > 1
  if (nterm > 1) {

    # create two parts, comma separated and `final_sep` separated
    first_part <- paste0(x[-nterm], collapse = ", ")
    second_part <- x[nterm]

    # return clean list
    x <- paste0(first_part, ", ", final_sep, " ", second_part)

  }

  # return clean list
  x

}

# update elements of an object based on a mask, function, and arguments
do_mask <- function(.x, mask, fun, ...) {
  args <- list(...)
  .x[mask] <- do.call(fun, c(list(.x[mask]), args))
  .x
}

# generate a unique ID for an object based on hex code
# function adapted from: https://github.com/greta-dev/greta
hex_id <- function() {
  paste(as.raw(sample.int(256L, 4, TRUE) - 1L), collapse = "")
}

# set an object class
as_class <- function(
  object,
  name,
  type = c("function", "list", "matrix", "array", "dynamics")
) {

  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))

  object

}

# check that a matrix is compatible with Leslie matrix updates
check_leslie <- function(x) {
  x <- x$matrix
  x[transition(x)] <- 0
  x[reproduction(x)] <- 0
  out <- TRUE
  if (any(x != 0))
    out <- FALSE
  out
}
