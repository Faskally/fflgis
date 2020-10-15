#' @export
getUpDownNodes <- function(rivs) {
  edges <- lapply(rivs @ lines, function(x) x @ Lines[[1]] @ coords)
  cc_up_nodes <- t(sapply(edges, function(x) head(x, 1)))
  cc_down_nodes <- t(sapply(edges, function(x) tail(x, 1)))

  cc2label <- function(x) paste0(x, collapse = ":")
  up_nodes <- apply(cc_up_nodes, 1, cc2label)
  down_nodes <- apply(cc_down_nodes, 1, cc2label)

  list(up_node = up_nodes, down_node = down_nodes)
}
