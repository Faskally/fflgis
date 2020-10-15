# Updated addSite2DRN function due to update to snapPointsToLines
#' @export
addSite2DRN <- function(site, rivs, site_name, rname,mdist = 200) {
  # find segment to split
  snap_site <- snapPointsToLines(site, rivs,maxDist = mdist)

  # is site already a node?
  if (paste(coordinates(snap_site), collapse = ":") %in% unlist(getUpDownNodes(rivs))) {
    # TODO  if node is network node add site data:
    message("skipping site ", paste(site_name, collapse = "; "))
    return(rivs)
  }

  # strip off river lines that we will keep
  # sometime there is more than one nearerst line id
  ids <- as.character(unlist(snap_site@data[,names(snap_site) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(snap_site, rivs[x,]))))
  }
  keep <- rownames(rivs@data) != ids
  rivs_keep <- rivs[keep,]
  # rename IDS
  spChFIDs(rivs_keep) <- 1:length(rivs_keep)

  # split segement:
  out <- c(cutLineDownstream(rivs[ids,], snap_site), cutLineUpstream(rivs[ids,], snap_site))
  # give unique ids and merge
  out <- lapply(1:2, function(i) {out[[i]]@lines[[1]]@ID <- paste(i); out[[i]]})
  out <- do.call(rbind, out)
  # change ids
  spChFIDs(out) <- 1:2 + length(rivs_keep)

  # get correct data.frame to go with new spatial lines
  newdf <- rivs[paste(snap_site$nearest_line_id),]@data[c(1,1),]
  if ("up_node" %in% names(newdf)) {
    newdf $ up_node <- getUpDownNodes(out)$up_node
    newdf $ down_node <- getUpDownNodes(out)$down_node
  }
  if ("start" %in% names(newdf)) {
    newdf $ start <- newdf $ up_node
    newdf $ end <- newdf $ down_node
  }
  newdf[[rname]][2] <- site_name

  row.names(newdf) <- 1:2 + length(rivs_keep)

  out <- SpatialLinesDataFrame(out, newdf)

  # now bind these segments onto the river
  rbind(rivs_keep, out)
}
