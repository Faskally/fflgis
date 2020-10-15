#' @export
plotPath <- function(path, g) {
  vgraph <- induced.subgraph(g, path)
  E(vgraph)$width <- 3
  E(vgraph)$color <- "orange"
  E(vgraph)$arrow.size <- 0.1
  V(vgraph)$label <- NA
  V(vgraph)$size <- 0.1

  plot(vgraph, add = TRUE, rescale = FALSE)
}
