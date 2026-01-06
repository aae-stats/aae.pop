#' @name dynamics
#' @title Create and update population dynamics objects
#' @description Define population dynamics from a matrix and additional
#'   objects that determine covariate effects, density dependence, and
#'   forms of stochasticity.
NULL

#' @rdname dynamics
#'
#' @export
#'
#' @param matrix a matrix of vital rates specifying transitions between
#'   ages or stages. Specified in the format $n_{t+1} = A %*% n_t$, where
#'   $A$ is the matrix, so that values in a given column and row denote a
#'   transition from that column to that row
#' @param \dots additional objects used to define population dynamics.
#'   Must be one or more of \code{\link{covariates}},
#'   \code{\link{replicated_covariates}},
#'   \code{\link{environmental_stochasticity}},
#'   \code{\link{demographic_stochasticity}},
#'   \code{\link{density_dependence}},
#'   \code{\link{add_remove_pre}}, or
#'   \code{\link{add_remove_post}}.
#'   Note that \code{\link{density_dependence_n}} is equivalent to
#'   \code{\link{add_remove_post}}.
#'
#' @details A call to \code{dynamics} defines an object of class
#'   \code{dynamics}, which can be used to simulate population
#'   trajectories with the \code{\link{simulate}} function. The \code{plot}
#'   function is supported and will generate a general life-cycle
#'   diagram based on the defined population dynamics.
#'
#' @examples
#' # define a population
#' nclass <- 5
#' popmat <- matrix(0, nrow = nclass, ncol = nclass)
#' popmat[reproduction(popmat, dims = 4:5)] <- c(10, 20)
#' popmat[transition(popmat)] <- c(0.25, 0.3, 0.5, 0.65)
#'
#' # define a dynamics object
#' dyn <- dynamics(popmat)
#'
#' # and plot this
#' plot(dyn)
dynamics <- function(matrix, ...) {
  # check matrix is provided
  if (missing(matrix)) {
    stop("matrix must be provided when creating dynamics object")
  }

  # is matrix actually a matrix?
  if (length(dim(matrix)) != 2) {
    stop("matrix must be a two-dimensional array or matrix ",
         "defining a population dynamics model",
         call. = FALSE
    )
  }

  # list possible processes
  processes_supported <- c(
    "covariates",
    "replicated_covariates",
    "environmental_stochasticity",
    "demographic_stochasticity",
    "density_dependence",
    "add_remove_pre",
    "add_remove_post",
    "density_dependence_n"
  )

  # collapse processes into a list
  processes <- list(...)

  # remove any named NULL processes
  processes <- processes[!sapply(processes, is.null)]

  # and work out which have been supplied
  processes_supplied <- sapply(processes, function(x) class(x)[1])
  names(processes) <- processes_supplied

  # error if processes not supported
  if (!all(processes_supplied %in% processes_supported)) {
    stop("Additional arguments to dynamics must be one of ",
         clean_paste(processes_supported, final_sep = "or"),
         call. = FALSE
    )
  }

  # error if multiple of the same process supplied
  if (length(processes_supplied) > 0) {
    nproc <- table(processes_supplied)
    if (any(nproc > 1)) {
      duplicate_process <- names(nproc)[nproc > 1]
      stop("Multiple objects provided for the following processes: ",
           clean_paste(duplicate_process, final_sep = "and"), ".\n",
           " A dynamics object can include up to one of each process type",
           call. = FALSE
      )
    }
  }

  # make a list of supplied processes, NULL if missing
  process_list <- processes[processes_supported]
  names(process_list) <- processes_supported

  # compile and return everything
  as_dynamics(
    c(
      list(
        nclass = nrow(matrix),
        hex = hex_id(),
        matrix = matrix
      ),
      process_list
    )
  )
}

#' @rdname dynamics
#'
#' @export
#'
#' @param object a \code{dynamics} object
#'
#' @details A compiled \code{dynamics} object can be
#'   updated to change any of the included processes with
#'   the \code{update} function.
update.dynamics <- function(object, ...) {
  # collate dots into a list
  processes_updated <- list(...)

  # remove any named NULL processes
  processes_updated <-
    processes_updated[!sapply(processes_updated, is.null)]

  # and work out which have been supplied
  processes_supplied <-
    sapply(processes_updated, function(x) class(x)[1])
  names(processes_updated) <- processes_supplied

  # pull out existing processes
  processes_supported <- c(
    "covariates",
    "replicated_covariates",
    "environmental_stochasticity",
    "demographic_stochasticity",
    "density_dependence",
    "add_remove_pre",
    "add_remove_post",
    "density_dependence_n"
  )
  processes_existing <- lapply(processes_supported, function(x) object[[x]])
  names(processes_existing) <- processes_supported

  # remove NULL (missing) processes
  processes_existing <-
    processes_existing[!sapply(processes_existing, is.null)]

  # update with new processes
  processes_existing[processes_supplied] <-
    processes_updated

  # recreate and return dynamics object with new processes
  do.call(dynamics, c(list(object$matrix), processes_existing))
}

# S3 plot method
#' @export
plot.template <- function(x, y, ...) {
  plot(x$dynamics)
}

# S3 plot method
#' @export
plot.dynamics <- function(
    x, y, ..., labels = NULL, cycle_first = "reproductive"
) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("the DiagrammeR package must be installed to ",
         "plot dynamics objects",
         call. = FALSE
    )
  }

  # pull out the population matrix
  mat <- x$matrix

  # work out the matrix structure (binary)
  mat <- ifelse(mat > 0, 1, 0)

  # and transpose it because it's easier than changing everything
  mat <- t(mat)

  # pull out process type
  type <- ifelse(any(diag(mat) > 0), "Stage", "Age")

  # add tweak for age-structured with final collecting stage
  if (all(diag(mat)[-nrow(mat)] == 0)) {
    type <- "Age"
  }

  gr <- DiagrammeR::from_adj_matrix(mat, # nolint
                                    mode = "directed",
                                    use_diag = TRUE
  )

  # how many nodes?
  n_nodes <- nrow(gr$nodes_df)

  # edge types
  to <- gr$edges_df$to
  from <- gr$edges_df$from

  # change type back to stage if there are multiple "age x+" terms
  if (sum(to == from) > 1) {
    type <- "Stage"
  }

  # identify different node types
  node_type <- rep("pre_reprod", n_nodes)
  if (cycle_first == "reproductive") {
    node_type[from[to == 1]] <- "reprod"
  } else {
    node_type[from[to == 1 & from > 1]] <- "reprod"
  }
  max_reprod <- max(which(node_type == "reprod"))
  node_type[seq_len(n_nodes) > max_reprod] <- "post_reprod"

  # change node types and colours based on node type
  node_shapes <- rep("circle", n_nodes)
  node_shapes[node_type == "reprod"] <- "circle"
  node_shapes[node_type == "post_reprod"] <- "circle"

  col_pal <- c("#57A0D3", "#4F7942", "#0E4D92", "#964000")
  col_pal_light <- paste0(col_pal, "50")
  node_edge_colours <- rep(col_pal[2], n_nodes)
  node_edge_colours[node_type == "reprod"] <- col_pal[3]
  node_edge_colours[node_type == "post_reprod"] <- col_pal[4]

  node_colours <- rep(col_pal_light[2], n_nodes)
  node_colours[node_type == "reprod"] <- col_pal_light[3]
  node_colours[node_type == "post_reprod"] <- col_pal_light[4]

  node_size <- rep(0.9, n_nodes)
  node_size[node_type == "reprod"] <- 0.9
  node_size[node_type == "post_reprod"] <- 0.9

  # add some labels for the nodes (age or stage depending on type of model),
  #    or user labels if provided
  if (is.null(labels)) {
    node_labels <- paste(type, seq_len(n_nodes), sep = " ")
  } else {

    # check that there are enough labels
    if (n_nodes != length(labels)) {
      stop(
        "There must be one label for each row/column in the population matrix",
        call = FALSE
      )
    }

    # use these if all OK
    node_labels <- labels

  }

  # if it's a Leslie matrix and to == from, we have an "age+" situation
  if (type == "Age" & is.null(labels)) {
    node_labels[from[to == from]] <- paste0(node_labels[from[to == from]], "+")
  }

  edge_style <- rep("solid", length(to))
  if (cycle_first == "reproductive") {
    edge_style[from %in% which(node_type == "reprod") & to == 1] <- "dashed"
  } else {
    edge_style[from %in% which(node_type == "reprod") &
                 !(to %in% which(node_type == "reprod")) &
                 !(to %in% which(node_type == "post_reprod"))] <- "dashed"
  }

  font_colour <- rep(col_pal[2], n_nodes)
  font_colour[node_type == "reprod"] <- col_pal[3]
  font_colour[node_type == "post_reprod"] <- col_pal[4]

  # node options
  gr$nodes_df$type <- node_type
  gr$nodes_df$fontcolor <- font_colour
  gr$nodes_df$fontsize <- 12
  gr$nodes_df$penwidth <- 2

  gr$nodes_df$shape <- node_shapes
  gr$nodes_df$color <- node_edge_colours
  gr$nodes_df$fillcolor <- node_colours
  gr$nodes_df$width <- node_size
  gr$nodes_df$height <- node_size * 0.8
  gr$nodes_df$label <- node_labels

  # edge options
  gr$edges_df$color <- "Gainsboro"
  gr$edges_df$fontname <- "Avenir"
  gr$edges_df$fontcolor <- "LightGray"
  gr$edges_df$fontsize <- 14
  gr$edges_df$penwidth <- 4

  edge_labels <- rep("transition", length(from))
  edge_labels <- ifelse(from == to, "survival", edge_labels)
  if (cycle_first == "reproductive") {
    edge_labels <- ifelse(from == to & from == 1, "reproduction", edge_labels)
  }
  edge_labels <- ifelse(from > to, "reproduction", edge_labels)

  gr$edges_df$label <- edge_labels
  gr$edges_df$style <- edge_style

  # set the layout type
  gr$global_attrs$value[gr$global_attrs$attr == "layout"] <- "dot"

  # make it horizontal
  gr$global_attrs <- rbind(
    gr$global_attrs,
    data.frame(
      attr = "rankdir",
      value = "LR",
      attr_type = "graph"
    )
  )

  gr_viz <- DiagrammeR::render_graph(gr) # nolint
  attr(gr_viz, "dgr_graph") <- gr
  gr_viz
}

# S3 is method
#' @rdname dynamics
#'
#' @export
#'
#' @param x an object to pass to \code{is.dynamics}
# nolint start
is.dynamics <- function(x) {
  # nolint end
  inherits(x, "dynamics")
}

# S3 print method
#' @export
# nolint start
print.dynamics <- function(x, ...) {
  # nolint end
  nsp <- x$nspecies
  if (is.null(nsp)) {
    nsp <- "a single"
  }

  cat(paste0("Population dynamics object for ", nsp, " species\n"))
}

# internal function: set dynamics class
as_dynamics <- function(x) {
  as_class(x, name = "dynamics", type = "list")
}
