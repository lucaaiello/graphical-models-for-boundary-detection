rm(list = ls())

# Standalone generation of main-manuscript Figure 4
# =============================================================================
# Run from the repository root with:
#   source("data analysis/disease_graph_schematics.R")
#
# This is a deterministic conceptual figure and does not use SEER observations
# or posterior draws. It reproduces the three four-disease graph topologies in
# the manuscript: unstructured (left), directed acyclic (center), and
# undirected (right). The output is:
#   data analysis/exploratory_figures/
#     figure_4_disease_graph_schematics.pdf

if (!all(file.exists(c("README.md", "data analysis")))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

output_dir <- file.path("data analysis", "exploratory_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_4_file <- file.path(
  output_dir, "figure_4_disease_graph_schematics.pdf"
)

node_coordinates <- rbind(
  `1` = c(0.25, 0.75),
  `2` = c(0.75, 0.75),
  `3` = c(0.25, 0.25),
  `4` = c(0.75, 0.25)
)
node_radius <- 0.095

complete_edges <- t(utils::combn(seq_len(nrow(node_coordinates)), 2L))
directed_edges <- rbind(
  c(1L, 2L), c(1L, 4L), c(2L, 3L), c(3L, 4L)
)
undirected_edges <- rbind(
  c(1L, 2L), c(1L, 4L), c(2L, 3L), c(3L, 4L)
)

edge_endpoints <- function(from, to, radius = node_radius) {
  start <- node_coordinates[from, ]
  end <- node_coordinates[to, ]
  direction <- end - start
  unit_direction <- direction / sqrt(sum(direction^2))
  list(
    start = start + radius * unit_direction,
    end = end - radius * unit_direction
  )
}

draw_graph_panel <- function(edges, directed = FALSE) {
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1), ylim = c(0, 1), asp = 1, xaxs = "i", yaxs = "i"
  )
  for (edge in seq_len(nrow(edges))) {
    endpoints <- edge_endpoints(edges[edge, 1L], edges[edge, 2L])
    if (directed) {
      graphics::arrows(
        endpoints$start[[1L]], endpoints$start[[2L]],
        endpoints$end[[1L]], endpoints$end[[2L]],
        length = 0.12, angle = 24, lwd = 1.4
      )
    } else {
      graphics::segments(
        endpoints$start[[1L]], endpoints$start[[2L]],
        endpoints$end[[1L]], endpoints$end[[2L]], lwd = 1.4
      )
    }
  }
  graphics::symbols(
    node_coordinates[, 1L], node_coordinates[, 2L],
    circles = rep(node_radius, nrow(node_coordinates)), inches = FALSE,
    add = TRUE, bg = "grey80", fg = "black", lwd = 1.5
  )
  graphics::text(
    node_coordinates[, 1L], node_coordinates[, 2L],
    labels = rownames(node_coordinates), cex = 1.4
  )
}

# RESULT: Main Figure 4 - disease-graph schematics -----------------------------
grDevices::cairo_pdf(figure_4_file, width = 12, height = 4.1, onefile = FALSE)
graphics::par(mfrow = c(1L, 3L), mar = c(0.2, 0.2, 0.2, 0.2))
draw_graph_panel(complete_edges)
draw_graph_panel(directed_edges, directed = TRUE)
draw_graph_panel(undirected_edges)
grDevices::dev.off()

if (!file.exists(figure_4_file) || file.info(figure_4_file)$size <= 0) {
  stop("Main Figure 4 was not created.", call. = FALSE)
}

cat(
  "Main Figure 4 created:\n  - ",
  normalizePath(figure_4_file, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)
