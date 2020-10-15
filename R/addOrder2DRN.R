#' @export
addOrder2DRN <- function(rivs, graph) {
  rivs@data$ order <- E(graph)$order
  rivs@data$ width <- E(graph)$width
  rivs@data$ color <- E(graph)$color
  # rivs@data $ start <- ends(graph, E(graph))[,1]
  # rivs@data $ end <- ends(graph, E(graph))[,2]
  rivs <- addUpDownNodes(rivs)
  rivs
}
