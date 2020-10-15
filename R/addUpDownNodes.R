#' @export
addUpDownNodes <- function(rivs) {
  upDownNodes <- getUpDownNodes(rivs)
  rivs$start <- upDownNodes$up_node
  rivs$end <- upDownNodes$down_node
  rivs$up_node <- upDownNodes$up_node
  rivs$down_node <- upDownNodes$down_node
  rivs
}
