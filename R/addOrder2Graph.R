#' @export
addOrder2Graph <- function(graph, order) {
  V(graph)$order <- order
  V(graph)$color <- V(graph)$order

  # assign orders to edges
  for (i in 1:max(V(graph)$order)) {
    if (sum(V(graph)$order == i) > 1) {
      E(graph)[from(V(graph)[V(graph)$order == i])]$order <- i
    }
  }

  E(graph)$color <- rev(terrain.colors(max(E(graph)$order)))[E(graph)$order]
  E(graph)$width <- E(graph)$order / max(E(graph)$order) * 7

  graph
}
