# create a neat list of terms for messages, warnings, and errors
clean_paste <- function(x, final_sep = "and") {

  # want number of terms so we can separate the final term
  nterm <- length(x)

  # create two parts, comma separated and `final_sep` separated
  first_part <- paste0(x[-nterm], collapse = ", ")
  second_part <- x[nterm]

  # return clean list
  paste0(first_part, ", ", final_sep, " ", second_part)

}

# set an object class
as_class <- function(
  object,
  name,
  type = c("function", "list", "matrix", "array")
) {

  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))

  object

}
