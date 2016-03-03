





#' @export
findBuffer <- function(p, up_distance = 100, width = 25, search_buffer = 200,
                       search_width = NULL, debug = FALSE)
{
  # expects to have the following objects available:
  #   rivs
  #   g
  #   wareas
  #   wlines
  #   dtm

  # setup
  if (is.null(search_width)) search_width <- c(7, 10, 15, 30, 50)

  # this is the main algorithm for finding buffers
  message("finding xy shift... ")
  xyshift <- findShift(p, search_buffer = search_buffer, rivs_buffer = 100, debug = debug)

  # subset rivs and shift
  # set bbox cover upstream distance*1.5 to be on the safe side
  bbox <- gBuffer(p, width = max(search_buffer, up_distance*1.5))
  if (debug) plot(bbox, border = grey(0.7))

  wk_rivs <- rivs[as.vector(rgeos::gIntersects(rivs, bbox, byid = TRUE)),]
  wk_rivs <- shift(wk_rivs, xyshift[1], xyshift[2])
  if (debug) lines(wk_rivs)

  # set sepa line buffer width based on river order
  # snap point to shifted river segment
  p_snap <- maptools::snapPointsToLines(p, wk_rivs) # no limit
  if (debug) points(p_snap, pch = 16, cex = 0.5, col = "red")

  # find search width first by strahler order:
  #   find line on river to get strahler order from
  ids <- as.character(unlist(p_snap@data[,names(p_snap) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap, wk_rivs[x,]))))
  }
  strahler <- wk_rivs[ids, ]$order
  search_width <- search_width[pmin(strahler, length(search_width))]

  message("finding buffer... ")
  # now get buffer

  # crop main datasets to size
  # use buffer to crop water lines
  wk_wlines <- wlines[as.vector(rgeos::gIntersects(wlines, bbox, byid = TRUE)),]
  wk_wlines <- gIntersection(wk_wlines, bbox)
  if (debug && !is.null(wk_wlines)) lines(wk_wlines, col = "blue")

  # use buffer to crop water areas
  wk_wareas <- wlines[as.vector(rgeos::gIntersects(wareas, bbox, byid = TRUE)),]
  wk_wareas <- gIntersection(wk_wareas, bbox)
  if (debug && !is.null(wk_wareas)) plot(wk_wareas, col = "blue", add = TRUE)

  wk_dtm <- crop(dtm, bbox)

  ## walk up sepa river
  out <- walkUpstream(p_snap, wk_rivs, up_distance)
  p_upstr <- out$p_upstr
  seg3 <- out$seg
  if (debug) {
    points(p_upstr, pch = 16, cex = 0.5, col = "orange")
    lines(seg3, lwd = 2)
  }
  # was p_snap a source?
  if (SpatialLinesLengths(seg3) == 0) {
    stop("The sample point snapped exactly to a source...")
  }

  # cut water polygon boundary here
  # tricky - one solution is
  buff_sepa <- rgeos::gBuffer(seg3, width = search_width, byid = FALSE, capStyle = "FLAT")
  if (debug) plot(buff_sepa, add = TRUE, border = grey(0.7))
  # start new plot
  if (debug) {
    plot(buff_sepa, border = grey(0.7), main = paste0("buffers using strahler: ", strahler))
    points(p_snap, col = "red", pch = 16, cex = 0.5)
    points(p_upstr, col = "orange", pch = 16, cex = 0.5)
    lines(seg3, col = grey(0.7))
  }

  # cut out water polygons MM data
  if (!is.null(wk_wareas)) {
    cut_wareas <- rgeos::gIntersection(buff_sepa, wk_wareas, byid=TRUE, drop_lower_td=TRUE)
    if (debug) plot(cut_wareas, col = "lightblue", add = TRUE)
  } else {
    cut_wareas <- NULL
  }

  # cut out water lines MM data
  if (!is.null(wk_wlines)) {
    cut_wlines <- rgeos::gIntersection(buff_sepa, wk_wlines, byid=TRUE, drop_lower_td=TRUE)
    if (debug) lines(cut_wlines, col = "blue")
  } else {
    cut_wlines <- NULL
  }

  if (is.null(cut_wareas) & is.null(cut_wlines)) {
    # there are no MM_water bodies near by... this needs to be flagged somehow...
    stop("no MM data for this SEPA river.")
  }

  # now add a buffer to the river segment
  buff_riva <- if (is.null(cut_wareas)) NULL else rgeos::gBuffer(cut_wareas, width = width, byid = FALSE)
  buff_rivl <- rgeos::gBuffer(cut_wlines, width = width, byid = FALSE)
  # join the buffers incase there is a slight discrepancy
  buff_riv <- if (is.null(buff_riva)) buff_rivl else  gUnion(buff_riva, buff_rivl)
  if (debug) plot(buff_riv, col = col_alpha("orange", 0.2), add = TRUE)

  # --------------------------------------------------------
  # --------------------------------------------------------
  #                       TODO
  # now we need to take this buffer and recalculate the river segment and wlines, warea etc
  #
  # --------------------------------------------------------
  # --------------------------------------------------------

  # substract off river areas
  buff_land <- if (is.null(wk_wareas)) buff_riv else rgeos::gDifference(buff_riv, wk_wareas)
  if (debug) plot(buff_land, col = col_alpha("lightgreen", 0.2), add = TRUE)

  # prepare outputs
  out <- list(buffer = buff_riv, buffer_nowater = buff_land,
              p_upstr = p_upstr, riv_seg = seg3,
              cut_area = cut_wareas, cut_lines = cut_wlines,
              buffer_sepa = buff_sepa)

  # add row IDs
  if ("SpatialPointsDataFrame" %in% is(p)) {
    out$buffer <- gUnaryUnion(out$buffer)
    out$buffer@polygons[[1]]@ID <- rownames(p@data)

    out$buffer_nowater <- gUnaryUnion(out$buffer_nowater)
    out$buffer_nowater@polygons[[1]]@ID <- rownames(p@data)

    if (!is.null(cut_wareas)) {
      out$cut_area <- gUnaryUnion(out$cut_area)
      out$cut_area@polygons[[1]]@ID <- rownames(p@data)
    }

    out$buffer_sepa <- gUnaryUnion(out$buffer_sepa)
    out$buffer_sepa@polygons[[1]]@ID <- rownames(p@data)

    out$riv_seg <- gLineMerge(out$riv_seg)
    out$riv_seg@lines[[1]]@ID <- rownames(p@data)

    if (!is.null(out$cut_wlines)) {
      mcut_lines <- try(gLineMerge(out$cut_lines), silent = TRUE)
      if(inherits(mcut_lines, "try-error")) {
        # weird error - not sure why this is happening...
        # ctm == 520 site == 16, siteID == 495
        out$cut_lines <- NULL
      } else {
        out$cut_lines <- mcut_lines
        out$cut_lines@lines[[1]]@ID <- rownames(p@data)
      }
    }
  }

  # return stuff
  out

}





#' @export
findShift <- function(p, search_buffer = 200, rivs_buffer = 100, debug = FALSE) {
  # find the xy offset that minimises the distance between the sepa river line
  # and the river bank / line feature
  bbox <- gBuffer(p, width = search_buffer)

  # find first draft of river lines to get a buffer
  wk_rivs <- rivs[as.vector(rgeos::gIntersects(rivs, bbox, byid = TRUE)),]
  buf <- buffer(wk_rivs, width = rivs_buffer)
  if (debug) plot(buf, col = col_alpha("orange", 0.2), main = "finding xyshift")

  # use buffer to crop rivs
  wk_rivs <- rivs[as.vector(rgeos::gIntersects(rivs, buf, byid = TRUE)),]
  wk_rivs <- gIntersection(wk_rivs, buf)
  if (debug) lines(wk_rivs)

  # use buffer to crop water lines
  wk_wlines <- wlines[as.vector(rgeos::gIntersects(wlines, buf, byid = TRUE)),]
  wk_wlines <- gIntersection(wk_wlines, buf)
  if (debug && !is.null(wk_wlines)) lines(wk_wlines, col = "blue")

  # get points to match with lines
  lxy <- spsample(wk_rivs, n = 100, type = "regular")
  if (debug) points(lxy, pch = 16, cex = 0.5, col = "red")

  # optimise
  xyshift <-
    if (!debug){
      optim(c(0,0), function(par, ...) {
        shifted_lxy <- shift(lxy, par[1], par[2])
        sum(rgeos::gDistance(shifted_lxy, wk_wlines, byid = TRUE))
      }, method = "BFGS")$par
    } else if (debug) {
      # if debug, then plot the shifted points at each iteration of the optimiser
      optim(c(0,0), function(par, ...) {
        shifted_lxy <- shift(lxy, par[1], par[2])
        points(shifted_lxy, pch = ".")
        sum(rgeos::gDistance(shifted_lxy, wk_wlines, byid = TRUE))
      }, method = "BFGS")$par
    }

  if (debug) {
    shifted_wk_rivs <- shift(wk_rivs, xyshift[1], xyshift[2])
    lines(shifted_wk_rivs, col = "purple", lwd = 2)
  }

  xyshift
}







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


