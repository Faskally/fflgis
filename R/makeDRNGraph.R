#' @export
makeDRNGraph <- function(rivs) {
  # Get edges from SEPA DRN and work out the nodes (top and tail the segments)
  updown_nodes <- getUpDownNodes(rivs)

  g <- make_graph(c(do.call(rbind, updown_nodes)), directed = TRUE)

  # add attributes to graph
  # location for nodes - this looks unnessisary, but its is safe.
  cc_vertices <- do.call(rbind, strsplit(V(g)$name, ":"))
  mode(cc_vertices) <- "numeric"
  V(g)$x <- cc_vertices[, 1]
  V(g)$y <- cc_vertices[, 2]

  # find mouths
  V(g)$color <- "blue"
  V(g)$color[degree(g, mode = "out") == 0] <- "red"

  # add lengths of edges
  edges <- lapply(rivs @ lines, function(x) x @ Lines[[1]] @ coords)
  E(g)$weight <- sapply(edges, LineLength, longlat = FALSE, sum = TRUE)
  g
}
