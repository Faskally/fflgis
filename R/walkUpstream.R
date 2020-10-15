#' @export
walkUpstream <- function(p_snap, wk_rivs, up_distance = 100, useRiverOrder = TRUE)
{
  # get segment that snapped point is on
  # make graph of river
  wk_g <- makeDRNGraph(wk_rivs)
  # add in vertex names to wk_rivs
  wk_rivs@data[c("start", "stop")] <- getUpDownNodes(wk_rivs)

  #
  ids <- as.character(unlist(p_snap@data[,names(p_snap) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap, wk_rivs[x,]))))
  }
  seg <- wk_rivs[ids, ]
  # cut the segment to start at the snapped point
  seg2 <- cutLineDownstream(seg, p_snap)
  # how long is the line
  seg_len <- SpatialLinesLengths(seg2) # 17.4, so we will encouter a split

  if (seg_len < up_distance) {
    # what are the next segments?
    upVs <- names(V(wk_g)[nei(as.character(seg@data$start), mode = "in")])
    # how many segments are there?
    if (length(upVs) == 0)
    {
      # we have hit a source - return a buffer along the available segement?
      warning("A source was encountered, and only a buffer around the available line was used.")
      p_upstr <- getPointUpstream(seg2, seg_len+0.01)
      seg3 <- seg2

    } else if (length(upVs) == 1)
    {
      ids <- get.edge.ids(wk_g, c(upVs[1], as.character(seg@data$start)))
      upseg <- wk_rivs[ids,]
      # what if we encounter another junction?
      p_upstr <- getPointUpstream(upseg, up_distance - SpatialLinesLengths(seg2))
      seg3 <- cutLineUpstream(upseg, p_upstr)
      # join lines
      seg3 <- gUnion(seg2, seg3)

    } else if (length(upVs) > 1)
    {
      ids <- get.edge.ids(wk_g, c(upVs[1], as.character(seg@data$start), upVs[2], as.character(seg@data$start)))
      upsegs <- wk_rivs[ids,]
      # is one a higher order? (only if useRiverOrder is TRUE)
      if (!useRiverOrder | upsegs$order[1] == upsegs$order[2]) {
        # walk up first
        p_upstr1 <- getPointUpstream(upsegs[1,], up_distance - SpatialLinesLengths(seg2))
        # cut line
        seg3.1 <- cutLineUpstream(upsegs[1,], p_upstr1)
        # walk up 2nd
        p_upstr2 <- getPointUpstream(upsegs[2,], up_distance - SpatialLinesLengths(seg2))
        # cut line
        seg3.2 <- cutLineUpstream(upsegs[2,], p_upstr2)
        # join lines
        # ## NOTE: this has some strange behaviour ... it strips of detail...
        seg3 <- gUnion(gUnion(seg2, seg3.1), seg3.2)
        # join points
        p_upstr <- gUnion(p_upstr1, p_upstr2)
      } else {
        # walk up mainstem
        upseg <- upsegs[which.max(upsegs$order),]
        # what if we encounter another junction?
        p_upstr <- getPointUpstream(upseg, up_distance - SpatialLinesLengths(seg2))
        seg3 <- cutLineUpstream(upseg, p_upstr)
        # join lines
        seg3 <- gUnion(seg2, seg3)
      }
    }
  } else {
    # No need to cross junctions
    # get a point upstream
    p_upstr <- getPointUpstream(seg2, up_distance)
    # cut line upstream at new point
    seg3 <- cutLineUpstream(seg2, p_upstr)
  }
  list(seg  = seg3, p_upstr = p_upstr)
}
