





#' @export
findBuffer <- function(p, up_distance = 100, width = 25, search_buffer = 200,
                       search_width = NULL, rivs_buffer = 100, debug = FALSE,
                       rivs, g, wareas, wlines, 
                       snap_to_shifted_sepaline = TRUE)
{
  # setup
  if (is.null(search_width)) search_width <- c(7, 10, 15, 30, 50)

  # this is the main algorithm for finding buffers
  message("finding xy shift... ")
  xyshift <- findShift(p, search_buffer = search_buffer, rivs_buffer = rivs_buffer, debug = debug,
                       rivs = rivs, g = g, wareas = wareas, wlines = wlines)

  # subset rivs and shift
  # set bbox cover upstream distance*1.5 to be on the safe side
  bbox <- gBuffer(p, width = max(search_buffer, up_distance*1.5))
  if (debug) plot(bbox, border = grey(0.7))

  wk_rivs_orig <- rivs[as.vector(rgeos::gIntersects(rivs, bbox, byid = TRUE)),]
  wk_rivs <- shift(wk_rivs_orig, xyshift[1], xyshift[2])
  if (debug) lines(wk_rivs)

  # snap point to (possibly) shifted river segment
  if (snap_to_shifted_sepaline) {
    # snap to shifted line
    p_snap <- snapPointsToLines(p, wk_rivs) # no limit
  } else {
    # snap to origional line and shift
    p_snap <- snapPointsToLines(p, wk_rivs_orig) # no limit
    p_snap <- shift(p_snap, xyshift[1], xyshift[2])
  }
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
  wk_wlines <- wlines[as.vector(gIntersects(wlines, bbox, byid = TRUE)),]
  if (length(wk_wlines) == 0) {
    wk_wlines <- NULL
  } else {
    wk_wlines <- gIntersection(wk_wlines, bbox)
  }
  if (debug && !is.null(wk_wlines)) lines(wk_wlines, col = "blue")

  # use buffer to crop water areas
  wk_wareas <- wareas[as.vector(gIntersects(wareas, bbox, byid = TRUE)),]
  if (length(wk_wareas) == 0) {
    wk_wareas <- NULL
  } else {
    wk_wareas <- gIntersection(wk_wareas, bbox)
  }
  if (debug && !is.null(wk_wareas)) plot(wk_wareas, col = "lightblue", add = TRUE)

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
    warning("The sample point snapped exactly to a source...")
  }

  # cut water polygon boundary here
  # tricky - one solution is
  buff_sepa <- gBuffer(seg3, width = search_width, byid = FALSE, capStyle = "FLAT")
  if (debug) plot(buff_sepa, add = TRUE, border = grey(0.7))
  # start new plot
  if (debug) {
    plot(buff_sepa, border = grey(0.7), main = paste0("buffers using strahler: ", strahler))
    scalebar(up_distance, divs = 4, type = "bar")
    points(p_snap, col = "red", pch = 16, cex = 0.5)
    points(p_upstr, col = "orange", pch = 16, cex = 0.5)
    lines(seg3, col = grey(0.7))
  }

  # cut out water polygons MM data
  if (!is.null(wk_wareas)) {
    cut_wareas <- gIntersection(buff_sepa, wk_wareas, byid=TRUE, drop_lower_td=TRUE)
    if (debug & !is.null(cut_wareas)) plot(cut_wareas, col = "lightblue", add = TRUE)
  } else {
    cut_wareas <- NULL
  }

  # cut out water lines MM data
  if (!is.null(wk_wlines)) {
    cut_wlines <- gIntersection(buff_sepa, wk_wlines, byid=TRUE, drop_lower_td=TRUE)
    if (debug & !is.null(cut_wlines)) lines(cut_wlines, col = "blue")
  } else {
    cut_wlines <- NULL
  }

  # now add a buffer to the river segment
  if (is.null(cut_wareas) & is.null(cut_wlines)) {
    # there are no MM_water bodies near by... this needs to be flagged somehow...
    buff_riv <- gBuffer(seg3, width = width, byid = FALSE)
  } else {
    buff_riva <- if (is.null(cut_wareas)) NULL else rgeos::gBuffer(cut_wareas, width = width, byid = FALSE)
    buff_rivl <- if (is.null(cut_wlines)) NULL else rgeos::gBuffer(cut_wlines, width = width, byid = FALSE)

    if (is.null(cut_wareas)) {
      buff_riv <- buff_rivl
    } else if (is.null(cut_wlines)) {
      buff_riv <- buff_riva
    } else {
      buff_riv <- gUnion(buff_riva, buff_rivl)
    }
  }

  if (debug) plot(buff_riv, col = col_alpha("orange", 0.2), add = TRUE)

  # --------------------------------------------------------
  # --------------------------------------------------------
  #                       TODO
  # now we need to take this buffer and recalculate the river segment and wlines, warea etc
  #
  # --------------------------------------------------------
  # --------------------------------------------------------

  # substract off river areas
  buff_land <- if (is.null(wk_wareas)) buff_riv else gDifference(buff_riv, wk_wareas)
  if (debug & !is.null(buff_land)) plot(buff_land, col = col_alpha("lightgreen", 0.2), add = TRUE)


  # prepare outputs
  out <- list(buffer = buff_riv, buffer_nowater = buff_land,
              p_upstr = p_upstr, riv_seg = seg3,
              cut_area = cut_wareas, cut_lines = cut_wlines,
              buffer_sepa = buff_sepa,
              site = p)

  # add row IDs
  if ("SpatialPointsDataFrame" %in% is(p)) {
    out$buffer <- gUnaryUnion(out$buffer)
    out$buffer@polygons[[1]]@ID <- rownames(p@data)

    if (!is.null(buff_land)) {
      out$buffer_nowater <- gUnaryUnion(out$buffer_nowater)
      out$buffer_nowater@polygons[[1]]@ID <- rownames(p@data)
    }

    if (!is.null(cut_wareas)) {
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

  # return stuff
  out

}





#' @export
findShift <- function(p, search_buffer = 200, rivs_buffer = 100, debug = FALSE,
                      rivs, g, wareas, wlines) {
  # find the xy offset that minimises the distance between the sepa river line
  # and the river bank / line feature
  bbox <- gBuffer(p, width = search_buffer)

  # find first draft of river lines to get a buffer
  wk_rivs <- rivs[as.vector(rgeos::gIntersects(rivs, bbox, byid = TRUE)),]
  buf <- gBuffer(wk_rivs, width = rivs_buffer)
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











# group buffers into spatial dataframes
#' @export
groupBufferList <- function(buffer_list) {

  get_spolydf <- function(what, buffer_list, sites) {
    outlist <- lapply(buffer_list, "[[", what)
    notnull <- which(!sapply(outlist, is.null))
    SpatialPolygonsDataFrame(do.call(rbind, outlist[notnull]),
                             sites@data[notnull,,drop=FALSE])
  }

  get_sldf <- function(what, buffer_list, sites) {
    outlist <- lapply(buffer_list, "[[", what)
    notnull <- which(!sapply(outlist, is.null))
    SpatialLinesDataFrame(do.call(rbind, outlist[notnull]),
                          sites@data[notnull,,drop=FALSE])
  }

  get_sptsdf <- function(what, buffer_list, sites) {
    outlist <- lapply(buffer_list, "[[", what)
    notnull <- which(!sapply(outlist, is.null))
    SpatialPointsDataFrame(do.call(rbind, outlist[notnull]),
                           sites@data[notnull,,drop=FALSE])
  }

  # DO IT!
  sites <- do.call(rbind, lapply(buffer_list, "[[", "site"))

  buffer <- get_spolydf("buffer", buffer_list, sites)
  buffer_nowater <- get_spolydf("buffer_nowater", buffer_list, sites)
  cut_area <- get_spolydf("cut_area", buffer_list, sites)

  cut_lines <- get_sldf("cut_lines", buffer_list, sites)
  riv_seg <- get_sldf("riv_seg", buffer_list, sites)

  list(buffer = buffer, buffer_nowater = buffer_nowater, cut_area = cut_area,
       cut_lines = cut_lines, riv_seg = riv_seg,
       sites = sites)
}

