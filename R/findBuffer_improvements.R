#' @export
findBuffer <- function(p, up_distance = 100, width = 25, sepaWidth = NULL, searchWidth = NULL, demo = FALSE)
{
  # this is the main algorithm for finding buffers
  message("finding xy shift... ")
  xyshift <- findShift(p, up_distance = 200, sepaWidth = sepaWidth)
  
  # calculate buffer using shifted sepa line
  #message("aligning sepa river to master map ... ")
  xy <- list(x = coordinates(p)[,1] + c(-200, 200),
             y = coordinates(p)[,2] + c(-200, 200))
  wk_rivs <- cropFeature(rivs, xy, buffer = 500)
  #wk_rivs_shift <- shift(wk_rivs, xyshift[1], xyshift[2])
    
  # set sepa line (the search) buffer width based on river order
  message("finding buffer... ")
  
  # need to find search width first
  # find closest line segment on shifted sepa network
  p_snap <- maptools::snapPointsToLines(p, wk_rivs) # no limit
  ids <- as.character(unlist(p_snap@data[,names(p_snap) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap, rivs[x,]))))
  }
  seg <- rivs[ids, ]
  strahler <- seg$order
  if (is.null(searchWidth)) {
    searchWidth <- c(7, 10, 15, 30, 50)[pmin(strahler, 5)]
  } else {
    searchWidth <- searchWidth[pmin(strahler, length(searchWidth))]
  }
  
  # now get buffer
  out <- getBuffer(p, up_distance = up_distance, width = width, sepaWidth = searchWidth, shift = xyshift, demo = demo)
  
  out
}


#' @export
findShift <- function(p, up_distance = 100, sepaWidth = NULL, useRiverOrder = TRUE) {
  # find the xy offset that minimises the distance between the sepa river line
  # and the river bank / line feature
  # need to first get a wide buffer ignoring river order
  # think about using down stream lines to help with alignment
  if (is.null(sepaWidth)) sepaWidth <- 100
  out <- getBuffer(p, up_distance = up_distance, width = 25, sepaWidth = sepaWidth, useRiverOrder = useRiverOrder)
  
  xy <- list(x = coordinates(p)[,1] + c(-200, 200),
             y = coordinates(p)[,2] + c(-200, 200))
  
  # get points to match with lines
  lxy <- spsample(out $ riv_seg, n = 20, type = "regular")
  
  # optimise
  opt <-
    optim(c(0,0), function(par, ...) {
      slxy <- shift(lxy, par[1], par[2])
      sxy <- snapPointsToLines(slxy, out$cut_lines)
      dist <- sum(sqrt(rowSums((coordinates(sxy) - coordinates(slxy))^2)))
      dist
    }, method = "BFGS")
  
  opt$par
}




#' @export
getBuffer <- function(p, up_distance = 100, width = 25, sepaWidth = 50, shift = c(0,0), useRiverOrder = TRUE, demo = FALSE) {
  # expects to have the following objects available:
  #   rivs
  #   g
  #   wareas
  #   wlines
  #   dtm
  
  # work within a +- 200m square about site location
  xy <- list(x = coordinates(p)[,1] + c(-200, 200),
             y = coordinates(p)[,2] + c(-200, 200))
  
  # crop main datasets to size
  wk_area <- cropFeature(wareas, xy, buffer = 500)
  wk_lines <- cropFeature(wlines, xy, buffer = 500)
  wk_rivs <- cropFeature(rivs, xy, buffer = 500)
  wk_dtm <- crop(dtm, extent(xy))
  
  # find closest point on sepa network
  p_snap <- maptools::snapPointsToLines(p, rivs) # no limit!
  
  ## walk up sepa river
  out <- walkUpstream(p_snap, up_distance, useRiverOrder)
  p_upstr <- out$p_upstr
  seg3 <- out$seg
  # was p_snap a source?
  if (SpatialLinesLengths(seg3) == 0) {
    stop("The sample point snapped exactly to a source...")
  }
  
  # cut water polygon boundary here
  # tricky - one solution is
  seg3_shift <- shift(gLineMerge(seg3), shift[1], shift[2])
  buff_sepa <- rgeos::gBuffer(seg3_shift, width = sepaWidth, byid = FALSE, capStyle = "FLAT")
  
  cut_area <- rgeos::gIntersection(buff_sepa, wk_area, byid=TRUE, drop_lower_td=TRUE)
  cut_lines <- rgeos::gIntersection(buff_sepa, wk_lines, byid=TRUE, drop_lower_td=TRUE)
  
  if (is.null(cut_area) & is.null(cut_lines)) {
    # there are no MM_water bodies near by... this needs to be flagged somehow...
    stop("no MM data for this SEPA river.")
  }
  
  # now add a buffer to the river segment
  buff_riva <- if (is.null(cut_area)) NULL else rgeos::gBuffer(cut_area, width = width, byid = FALSE)
  
  # get buffer based on OS lines
  buff_rivl <- rgeos::gBuffer(cut_lines, width = width, byid = FALSE)
  
  # join the buffers incase there is a slight discrepancy
  buff_riv <- if (is.null(buff_riva)) buff_rivl else  gUnion(buff_riva, buff_rivl)
  
  # substract off river areas
  buff_land <- gDifference(buff_riv, wk_area)
  
  # prepare outputs
  out <- list(buffer = buff_riv, buffer_nowater = buff_land,
              p_upstr = p_upstr, riv_seg = seg3,
              cut_area = cut_area, cut_lines = cut_lines,  # gLineMerge(seg3)?
              buffer_sepa = buff_sepa)
  
  # add row IDs
  if ("SpatialPointsDataFrame" %in% is(p)) {
    out$buffer <- gUnaryUnion(out$buffer)
    out$buffer@polygons[[1]]@ID <- rownames(p@data)
    
    out$buffer_nowater <- gUnaryUnion(out$buffer_nowater)
    out$buffer_nowater@polygons[[1]]@ID <- rownames(p@data)
    
    if (!is.null(cut_area)) {
      out$cut_area <- gUnaryUnion(out$cut_area)
      out$cut_area@polygons[[1]]@ID <- rownames(p@data)
    }
    
    out$buffer_sepa <- gUnaryUnion(out$buffer_sepa)
    out$buffer_sepa@polygons[[1]]@ID <- rownames(p@data)
    
    out$riv_seg <- gLineMerge(out$riv_seg)
    out$riv_seg@lines[[1]]@ID <- rownames(p@data)
    
    if (!is.null(out$cut_lines)) {
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
  
  
  # plot results
  if (demo) {
    plot(wk_dtm)
    plot(wk_area, add = TRUE, border = NA, col = col_alpha("lightblue", 0.5))
    plot(wk_rivs, add = TRUE, col = "red", lwd = 2)
    plot(wk_lines, add = TRUE, col = "blue")
    points(p_snap, col = "gold", pch = 16, cex = 2)
    points(p_upstr, col = "orange", pch = 16, cex = 2)
    plot(seg3, col = "brown", add = TRUE, lwd = 3)
    plot(buff_sepa, col = paste0(grey(0.5), "22"), add = TRUE)
    plot(cut_area, add = TRUE, col = "purple")
    plot(cut_lines, add = TRUE, col = "purple", lwd = 2)
    plot(buff_riva, add = TRUE, col = col_alpha(grey(0.5), 0.5) , border = "green", lwd = 2)
    plot(cut_area, col = "red", add = TRUE)
    plot(buff_rivl, add = TRUE, col = col_alpha(grey(0.5), 0.5), border = "green", lwd = 2)
    plot(cut_lines, col = "pink", add = TRUE)
    plot(buff_riv, add = TRUE, col = col_alpha(grey(0.5), 0.5), border = "cyan", lwd = 2)
  }
  
  # return stuff
  out
}


#' @export
walkUpstream <- function(p_snap, up_distance = 100, useRiverOrder = TRUE)
{
  # get segment that snapped point is on
  ## expects to have the following objects available:
  #  rivs
  #  g
  
  ids <- as.character(unlist(p_snap@data[,names(p_snap) == "nearest_line_id"]))
  if (length(ids) > 1) {
    ids <- names(which.min(sapply(ids, function(x) gDistance(p_snap, rivs[x,]))))
  }
  seg <- rivs[ids, ]
  # cut the segment to start at the snapped point
  seg2 <- cutLineDownstream(seg, p_snap)
  # how long is the line
  seg_len <- SpatialLinesLengths(seg2) # 17.4, so we will encouter a split
  if (seg_len < up_distance) {
    # what are the next segments?
    upVs <- names(V(g)[nei(as.character(seg@data$start), mode = "in")])
    # how many segments are there?
    if (length(upVs) == 0) {
      # we have hit a source - return a buffer along the available segement?
      warning("A source was encountered, and only a buffer around the available line was used.")
      p_upstr <- getPointUpstream(seg2, seg_len+0.01)
      seg3 <- seg2
    } else 
    if (length(upVs) == 1) {
      ids <- get.edge.ids(g, c(upVs[1], as.character(seg@data$start)))
      upseg <- rivs[ids,]
      # what if we encounter another junction?
      p_upstr <- getPointUpstream(upseg, up_distance - SpatialLinesLengths(seg2))
      seg3 <- cutLineUpstream(upseg, p_upstr)
      # join lines
      seg3 <- gUnion(seg2, seg3)
    } else {
      ids <- get.edge.ids(g, c(upVs[1], as.character(seg@data$start), upVs[2], as.character(seg@data$start)))
      upsegs <- rivs[ids,]
      # is one a higher order? (only if useRiverOrder is TRUE)
      if (!useRiverOrder | upsegs$order[1] == upsegs$order[2]) {
        # walk up first
        p_upstr1 <- getPointUpstream(upsegs[1,], up_distance - SpatialLinesLengths(seg2))
        # cut line
        seg3.1 <- cutLineUpstream(upsegs[2,], p_upstr1)
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
