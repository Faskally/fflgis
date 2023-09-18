#' @export
getSeaDistance <- function(p, rivs, g) {
  p_snap <- snapPointsToLines(p, rivs)

  ids <- as.character(unlist(p_snap@data[,names(p_snap) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap, rivs[x,]))))
  }
  seg <- rivs[ids, ]
  # cut the segment to start at the snapped point
  seg2 <- cutLineUpstream(seg, p_snap)
  # how long is the line
  seg_len <- SpatialLinesLengths(seg2)
  # now find the downstream node name
  # convert segment to 1 edged graph and find mouth
  gs <- makeDRNGraph(seg)
  downV <- V(g)[names(summariseDRN(gs)$mouth)]

  # get the shortest path to the mouth
  vmouth <- summariseDRN(g) $ mouth
  vpath <- get.shortest.paths(g, from = downV, to = vmouth) $ vpath[[1]]
  epath <- E(g)[vpath]

  # get total distance
  dist <- seg_len + sum(epath$weight)

  list(distance = dist, path = vpath)
}


#' @export
getSeaDistance1 <- function(p_snap, seg) {

  # cut the segment to start at the snapped point
  seg2 <- cutLineUpstream(seg, p_snap)
  # how long is the line
  seg_len <- SpatialLinesLengths(seg2)
  # now find the downstream node name
  # convert segment to 1 edged graph and find mouth
  updown_nodes <- getUpDownNodes(seg)

  list(downV = updown_nodes$down_node, seg_len = seg_len)
}

#' @export
getSeaDistance2 <- function(downV, g) {

  # get the shortest path to the mouth
  vmouth <- summariseDRN(g) $ mouth
  vpath <- get.shortest.paths(g, from = downV, to = vmouth) $ vpath[[1]]
  epath <- E(g)[vpath]

  # get total distance
  dist <- sum(epath$weight)

  dist
}
