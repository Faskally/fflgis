#' @export
stripOrder <- function(graph) {

  # set local options
  opt <- igraph_opt("return.vs.es")
  igraph.options(return.vs.es = TRUE)
  on.exit(igraph.options(return.vs.es = opt))

  stopifnot(is.named(graph))
  wk_graph <- graph

  # lop while there are psuedo source nodes next to sources
  while(length(psourceID <- idPsuedoSource(wk_graph, graph)) > 0) {
    #print(psourceID)
    wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% psourceID)
  }
  # All remaining outer nodes are to be tagged and pruned
  outerID <- names(which(degree(wk_graph, mode = "in") == 0))
  if (length(outerID)==0) {
    message("Possible problem with graph: no source nodes but nvertices > 1\n",
            "This can happen if there is a cycle in the graph, or an edge pointing backwards.")
    return(wk_graph)
  }
  wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% outerID)

  wk_graph
}
