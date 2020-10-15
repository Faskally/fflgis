#' @export
getOrder <- function(graph, plot.it = FALSE, max_order = 10) {

  # set local options
  opt <- igraph_opt("return.vs.es")
  igraph.options(return.vs.es = TRUE)
  on.exit(igraph.options(return.vs.es = opt))

  stopifnot(is.named(graph))
  wk_graph <- graph

  i <- 1
  order <- rep(NA, vcount(graph))
  names(order) <- V(graph)$name

  if (plot.it && vcount(wk_graph) > 0) {
    plot(wk_graph,
      vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1,
      main = i, rescale = FALSE,
      xlim = range(V(graph)$x), ylim = range(V(graph)$y)
    )
  }

  while (vcount(wk_graph) > 1) {
    if (i > max_order) {
      message(
        "Possible problem with graph: maximum order reached.  \n",
        "Returning graph."
      )
      return(wk_graph)
    }

    # print(i)
    # lop while there are psuedo source nodes next to sources
    while (length(psourceID <- idPsuedoSource(wk_graph, graph)) > 0) {
      # print(psourceID)
      order[psourceID] <- i
      wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% psourceID)
    }
    # All remaining outer nodes are to be tagged and pruned
    outerID <- names(which(degree(wk_graph, mode = "in") == 0))
    if (length(outerID) == 0) {
      message(
        "Possible problem with graph: no source nodes but nvertices > 1\n",
        "This can happen if there is a cycle in the graph, or an edge pointing backwards.\n",
        "Returning graph."
      )
      return(wk_graph)
    } # else if (length(outerID)==1) {
    # then we are on the max order?
    # }
    order[outerID] <- i
    wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% outerID)

    if (plot.it && vcount(wk_graph) > 0) {
      plot(wk_graph,
        vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1,
        main = i + 1, rescale = FALSE,
        xlim = range(V(graph)$x), ylim = range(V(graph)$y)
      )
    }

    i <- i + 1
  }
  if (sum(is.na(order)) == 1) order[is.na(order)] <- i
  ### DONE

  order
}
