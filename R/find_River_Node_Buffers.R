# Edited buffer function for river nodes

find_River_Node_Buffers<-function (p, up_distance = 100, width = 25, search_buffer = 200,
                                   search_width = NULL, rivs_buffer = 100, debug = FALSE, rivs,
                                   g, wareas, wlines, snap_to_shifted_sepaline = TRUE)
{
  #browser()
  # This function is slightly different from Colins findBuffers function because the
  # findBuffers function was originally designed around finding buffers for sites which
  # we added into the network. Yet when you are finding buffers for nodes in the network, they
  # are already snapped to the lines and thus the original findBuffers function was suboptimal because:
  # 1. it was not reproduceable because the xyshift was different each time, the nearest_line_id i.e.
  # the first line segment in the distance matrix closest to the node could be different each time
  # and thus the buffers different each time - also these buffers might not always be correct
  # 2. ideally we would want to be able to specify that rather than just using which.min within
  # snapPointsToLines the closest line segment should be the segment with the highest order / or
  # nearest coordinates on the end node so that the buffers would always start on a mainstem and
  # not a tributary. For example if you had two tributaries and the nearest_line_ID from your
  # snapPointsToLines was one of the tributaries the code would never know about the other tributary when
  # it looked upstream. Alternatively if it snapped to the end node of the mainsteam/highest order it would
  # "see" two line segments going up from it.
  # These two issues are explained further and should become clearer in the subsequent annotation of this function :-)
  if (is.null(search_width))
    search_width <- c(7, 10, 15, 30, 50)
  message("finding xy shift... ")
  # This is a function which looks to move the drn the smallest amount to create the best match
  # with the MM data - we do this to ensure we can clip out as much of the wareas and wlines despite
  # the differences in digitisation
  # HOWEVER - this xyshift is different each time because it uses "optim" within the function, which is
  # an optimiser for finding the smallest shift. Thus this could cause sites/rivers to shift to
  # slightly different locations each time which will have implications for the buffers when these
  # nodes are already on confluences i.e. prediction nodes (this would not happen for sites added to the network)
  xyshift <- findShift(p, search_buffer = search_buffer, rivs_buffer = rivs_buffer,
                       debug = debug, rivs = rivs, g = g, wareas = wareas, wlines = wlines)
  # print(xyshift)
  # bbox is just a radius around the site to allow you to clip the datasets
  bbox <- gBuffer(p, width = max(search_buffer, up_distance *
                                   1.5))
  if (debug)
    plot(bbox, border = grey(0.7))
  # Clipped original drn
  wk_rivs_orig <- rivs[as.vector(rgeos::gIntersects(rivs, bbox,
                                                    byid = TRUE)), ]
  #print(wk_rivs_orig@data)
  # Shifted drn (using xyshift so can be slightly different each time)
  wk_rivs <- shift(wk_rivs_orig, xyshift[1], xyshift[2])
  #print(wk_rivs@data)
  if (debug)
    lines(wk_rivs)
  if (snap_to_shifted_sepaline) {
    # These section is completely different to to findBuffers code which did : p_snap <- snapPointsToLines(p, wk_rivs)
    # We do not want to snap the points to wk_rivs (i.e. the shifted river) though because it is slightly different
    # each time which is why we may get different buffers for river nodes. Instead we do this:
    # Snap to the original non shifted river but first do 0 shift on it - the only reason we are doing it
    # is so that the nearest_line_ID from snapPointsToLines then becomes the row position in the drn
    # rather than the rowname if you do not do the shift - snapPointsToLines will return the
    # row number (e.g. "127") and then the code will fall over later at gDistance(p_snap,
    # wk_rivs[x, ] when looks for say row number x "127" when it really wants "2" - the second line on drn@data
    # Thus the fake shift means we have the line position not the name
    wk_rivs_orig <- shift(wk_rivs_orig, 0, 0)
    # Find the line which is closes to the point - now for added sites choosing the minimum distance is not a problem
    # However for nodes which are already part of the network this is problematic because you could have more than one
    # river segment which is 0 from the point i.e. one river line where the node coordinates are the upstream point (mainstem)
    # and then two where the node coordinates are the downstream points (tributaries). This is thus problematic when which.min
    # is used because which.min just chooses the FIRST smallest distance in the distances matrix - which could be different
    # each time depending on the order of the matrix (which could change when the xyshift is different and the node snaps to
    # different nearest lines)
    p_snap <- snapPointsToLinesbyOrder(p, wk_rivs_orig)
    # to understand the snapPointsTOLines function (from maptools) and how we would like to edit it - here is a small description
    # of the process used within it:
    # 1. get a matrix of distances of all the points to all this lines within wk_rivs
    # 2. get a vector of the "nearest_line_index" which is the row in the dataframe which contains the line
    # segment nearest to the node/site - IMPORTANTLY which.min chooses the first in the list which means that where there
    # are more than one line segments which are 0 from the node (as explained above - row 59 onwards) only the first is chosen
    # this could be any of the lines - although by using wk_rivs_org we ensure it is the same line each time - whereas when snapping
    # to the shifted lines with shifted points this may not be the case
    # 3. adds the column nearest_line_id to the points@data which is the ID of the line which is closest to the node
    # row of the drn@data which corresponds to the ID of the lines bit of the sldf (drn@lines$ID)

    # HOW WE WANT TO EDIT IT - Instead of picking the first instance we want to pick the river segment with the highest order /
    # where the coordinates match the upnode because this will ensure we are always starting the buffer from the mainstem rather
    # than a tributary which will ensure that the buffer always goes up both tributaries if they are of the same order
    # SEE ANNOTATION IN THE snapPointsToLinesbyOrder function
  }
  else {
    # Do not want to use the below processing on the river nodes because the xyshift is slightly difference each time
    # which means the nodes may snap to different river segments each time - originally all 0 distances away.
    # Snap to original river
    p_snap <- snapPointsToLines(p, wk_rivs_orig)
    #print(p_snap@data)
    # shift the point/node by the xyshift (this xyshift is slightly different each time so could move to
    # a different river segment each time)
    p_snap <- shift(p_snap, xyshift[1], xyshift[2])
    #print(p_snap@data)
    # remove the original nearest line id
    p_snap <- p_snap[, !names(p_snap) %in% "nearest_line_id"]
    #print(p_snap@data)
    # snap the shift points/nodes to the shifted rivers (could be slightly different each time due to differences in xyshift
    # this would only be a problem where node is at a confluence and could snap to multiple lines - sites added into the
    # network should not have this problem)
    p_snap <- snapPointsToLines(p_snap, wk_rivs)
    #print(p_snap@data)
  }
  if (debug)
    points(p_snap, pch = 16, cex = 0.5, col = "red")
  # get the ids (which are row positions) for the nearest line to each node/site
  ids <- as.character(unlist(p_snap@data[, names(p_snap) ==
                                           "nearest_line_id"]))
  if (length(ids) > 1) {
    # if there is more than one - want to minimum distance - again which.min just picks the first instance
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap,
                                                             wk_rivs[x, ]))))
  }
  # get the river order for each nearest line
  strahler <- wk_rivs[ids, ]$order
  # search widths for each river line - the search width is different depending on orderand specified in the
  # function call
  search_width <- search_width[pmin(strahler, length(search_width))]
  message("finding buffer... ")
  # extract water lines around the node
  wk_wlines <- wlines[as.vector(gIntersects(wlines, bbox, byid = TRUE)),
                      ]
  if (length(wk_wlines) == 0) {
    wk_wlines <- NULL
  }
  else {
    wk_wlines <- gIntersection(wk_wlines, bbox)
  }
  if (debug && !is.null(wk_wlines))
    lines(wk_wlines, col = "blue")
  # extract water areas
  wk_wareas <- wareas[as.vector(gIntersects(wareas, bbox, byid = TRUE)),
                      ]
  if (length(wk_wareas) == 0) {
    wk_wareas <- NULL
  }
  else {
    wk_wareas <- gIntersection(wk_wareas, bbox)
  }
  if (debug && !is.null(wk_wareas))
    plot(wk_wareas, col = "lightblue", add = TRUE)
  # walk upstream the specified distance of the point
  # WalkUpstream is another function:
  out <- walkUpstream(p_snap, wk_rivs, up_distance)
  # Get upstream point
  p_upstr <- out$p_upstr
  # Get river /drn segment in the upstream distance
  seg3 <- out$seg
  if (debug) {
    points(p_upstr, pch = 16, cex = 0.5, col = "orange")
    lines(seg3, lwd = 2)
  }
  if (SpatialLinesLengths(seg3) == 0) {
    warning("The sample point snapped exactly to a source...")
  }
  # Buffer around the river segment
  buff_sepa <- gBuffer(seg3, width = search_width, byid = FALSE,
                       capStyle = "FLAT")
  if (debug)
    plot(buff_sepa, add = TRUE, border = grey(0.7))
  if (debug) {
    plot(buff_sepa, border = grey(0.7), main = paste0("buffers using strahler: ",
                                                      strahler))
    scalebar(up_distance, divs = 4, type = "bar")
    points(p_snap, col = "red", pch = 16, cex = 0.5)
    points(p_upstr, col = "orange", pch = 16, cex = 0.5)
    lines(seg3, col = grey(0.7))
  }
  if (!is.null(wk_wareas)) {
    # Buffer around the water areas (this can be null if the MM location is only digitised by lines)
    cut_wareas <- gIntersection(buff_sepa, wk_wareas, byid = TRUE,
                                drop_lower_td = TRUE)
    if (debug & !is.null(cut_wareas))
      plot(cut_wareas, col = "lightblue", add = TRUE)
  }
  else {
    cut_wareas <- NULL
  }
  if (!is.null(wk_wlines)) {
    # Buffer around the water lines
    cut_wlines <- gIntersection(buff_sepa, wk_wlines, byid = TRUE,
                                drop_lower_td = TRUE)
    if (debug & !is.null(cut_wlines))
      lines(cut_wlines, col = "blue")
  }
  else {
    cut_wlines <- NULL
  }
  if (is.null(cut_wareas) & is.null(cut_wlines)) {
    # Buffer around the drn
    buff_riv <- gBuffer(seg3, width = width, byid = FALSE)
  }
  else {
    buff_riva <- if (is.null(cut_wareas))
      NULL
    else rgeos::gBuffer(cut_wareas, width = width, byid = FALSE)
    buff_rivl <- if (is.null(cut_wlines))
      NULL
    else rgeos::gBuffer(cut_wlines, width = width, byid = FALSE)
    if (is.null(cut_wareas)) {
      buff_riv <- buff_rivl
    }
    else if (is.null(cut_wlines)) {
      buff_riv <- buff_riva
    }
    else {
      buff_riv <- gUnion(buff_riva, buff_rivl)
    }
  }
  if (debug)
    plot(buff_riv, col = col_alpha("orange", 0.2), add = TRUE)
  buff_land <- if (is.null(wk_wareas))
    buff_riv
  else gDifference(buff_riv, wk_wareas)
  if (debug & !is.null(buff_land))
    plot(buff_land, col = col_alpha("lightgreen", 0.2), add = TRUE)
  # Build up everything that will be within the buffer list
  out <- list(buffer = buff_riv, buffer_nowater = buff_land,
              p_upstr = p_upstr, riv_seg = seg3, cut_area = cut_wareas,
              cut_lines = cut_wlines, buffer_sepa = buff_sepa, site = p)
  if ("SpatialPointsDataFrame" %in% is(p)) {
    # If there are multiple buffers (like two tributaries) merge them into one
    out$buffer <- gUnaryUnion(out$buffer)
    out$buffer@polygons[[1]]@ID <- rownames(p@data)
    if (!is.null(buff_land)) {
      # If there are multiple no water buffers (like two tributaries) merge them into one
      out$buffer_nowater <- gUnaryUnion(out$buffer_nowater)
      out$buffer_nowater@polygons[[1]]@ID <- rownames(p@data)
    }
    if (!is.null(cut_wareas)) {
      # If there are multiple cut areas (like two tributaries) merge them into one
      out$cut_area <- gUnaryUnion(out$cut_area)
      out$cut_area@polygons[[1]]@ID <- rownames(p@data)
    }
    # If there are multiple buffers (like two tributaries) merge them into one
    out$buffer_sepa <- gUnaryUnion(out$buffer_sepa)
    out$buffer_sepa@polygons[[1]]@ID <- rownames(p@data)
    out$riv_seg <- gLineMerge(out$riv_seg)
    out$riv_seg@lines[[1]]@ID <- rownames(p@data)
    if (!is.null(out$cut_lines)) {
      # If there are multiple cutlines (like two tributaries) merge them into one
      mcut_lines <- try(gLineMerge(out$cut_lines), silent = TRUE)
      if (inherits(mcut_lines, "try-error")) {
        out$cut_lines <- NULL
      }
      else {
        out$cut_lines <- mcut_lines
        out$cut_lines@lines[[1]]@ID <- rownames(p@data)
      }
    }
  }
  out
}
