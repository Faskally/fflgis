

#' @export
cutLineDownstream <- function(line, pt) {
  x <- line@lines[[1]]@Lines[[1]]
  #get digitised points along river
  cc = coordinates(x)
  # find the distance between the snapped point and the digitised points
  dists <- sqrt(colSums((t(cc) - c(coordinates(pt)))^2))
  # find which interval the snapped point is in
  int <- min(order(dists)[1:2])
  # new coords for line (cut downstream points out)
  cc_new <- rbind(cc[1:int,],coordinates(pt))
  SpatialLines(list(Lines(list(Line(cc_new)), ID = "A")), CRS(proj4string(line)))
}


#' @export
cutLineUpstream <- function(line, pt) {
  x <- line@lines[[1]]@Lines[[1]]
  #get digitised points along river
  cc = coordinates(x)
  # find the distance between the snapped point and the digitised points
  dists <- sqrt(colSums((t(cc) - c(coordinates(pt)))^2))
  # find which interval the snapped point is in
  int <- max(order(dists)[1:2])
  # new coords for line (cut upstream points out)
  cc_new <- rbind(coordinates(pt), cc[-(1:(int-1)),])
  SpatialLines(list(Lines(list(Line(cc_new)), ID = "A")), CRS(proj4string(line)))
}

#' @export
getPointUpstream <- function(line, dist) {
  x <- line@lines[[1]]@Lines[[1]]
  #get digitised points along river
  cc = coordinates(x)
  # get distances between digitised points
  lengths <- LineLength(cc, longlat = FALSE, sum = FALSE)
  # remove duplicate points on line
  if (any(abs(lengths) < .Machine$double.eps)) {
    wl <- which(abs(lengths) < .Machine$double.eps)
    cc <- cc[-(wl), ]
    lengths <- lengths[-(wl)]
  }
  # cumulative distance
  csl = c(0, cumsum(lengths))
  # max distance
  maxl = csl[length(csl)]
  # choose a distance upstream
  pts <- maxl - dist
  # find which subsegment pts is in
  int = findInterval(pts, csl, all.inside = TRUE)
  # where is pts on the subsegment
  where = (pts - csl[int])/diff(csl)[int]
  # work out coords for the new point
  xy = cc[int, , drop = FALSE] +
    where * (cc[int + 1, , drop = FALSE] - cc[int, , drop = FALSE])
  SpatialPoints(xy, CRS(proj4string(line)))
}


#' @export
cropFeature <- function(x, xy, buffer = 0) {
  # make a spatial polygon for the new limits
  bbox <- cbind(x = sort(xy$x)[c(1,2,2,1,1)],
                y = sort(xy$y)[c(1,1,2,2,1)])
  bbox <- SpatialPolygons(list(Polygons(list(Polygon(bbox)), ID = "a")))
  bbox @ proj4string <- x @ proj4string

  if (buffer > 0) {
    bbox <- gBuffer(bbox, width = buffer, byid = TRUE)
  }

  ind <- as.vector(rgeos::gIntersects(x, bbox, byid = TRUE))

  if (sum(ind) == 0) ind <- 1

  x[ind,]
}




#' @export
idPsuedoSource <- function(gr, gr_full) {
  # determine if the outermost nodes are to be pruned back
  sourceID <- character(0)
  # get all psuedo nodes
  psuedoID <- names(which(degree(gr, mode = "in") == 1))
  # identify if these nodes are one down from a source
  if (length(psuedoID) > 0) {
    # get upstream neighbours
    upID <- sapply(neighborhood(gr, 1, psuedoID, mode = "in"), function(x) names(x)[2])
    # is the upID a source, if yes mark it to be removed
    sourceID <- c(sourceID, unique(names(which(degree(gr, mode = "in")[upID] == 0))))
  }
  # get all braided ends
  confluenceID <- names(which(degree(gr, mode = "in") > 1))
  if (length(confluenceID) > 0) {
    # what are the upstream nodes
    upID <- lapply(neighborhood(gr, 1, confluenceID, mode = "in"), function(x) names(x)[-1])
    # are these upstream nodes sources?
    is.PsuedoSource <- sapply(upID, function(id) all(degree(gr, mode = "in")[id] == 0))
    # if so are these upstream nodes connected?
    if (any(is.PsuedoSource)) {
      is.braid <- sapply(upID[is.PsuedoSource], function(x) {
        if (length(x) == 1) return(TRUE) # odd one... why would this happen??
        # how big a neighbourhood do we consider for a braid to be a braid?
        # I have chosen 20 here... but it could easily be argued to be something less than 10...
        upIDs <- lapply(neighborhood(gr_full, 20, x, mode = "in"), names)
        if (length(x) == 2) {
          length(do.call(intersect, upIDs)) > 0
        } else {
          # if there is atleast one unconnected node then this is
          # a confulence,
          # otherwise if all up nodes are connected then it is a braid.
          upcommon <- apply(combn(1:length(upIDs[c(1,1,2)]), 2), 2, function(ids) length(do.call(intersect, upIDs[ids])))
          all(upcommon > 0)
        }
      })
      # mark confulenceIDs that are braids as sources to be removed
      sourceID <- c(sourceID, unlist(upID[is.PsuedoSource][is.braid]))
    }
  }

  sourceID
}


#' @export
getOrder <- function(graph, plot.it = FALSE) {

  stopifnot(is.named(graph))
  wk_graph <- graph

  i <- 1
  order <- rep(NA, vcount(graph))
  names(order) <- V(graph)$name

  if (plot.it && vcount(wk_graph) > 0) {
    plot(wk_graph, vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1, main = i)
  }

  while(vcount(wk_graph) > 1) {
    if (i > max_order) {
      message("Possible problem with graph: maximum order reached.")
      return(order)
    }
  
    # lop while there are psuedo source nodes next to sources
    while(length(psourceID <- idPsuedoSource(wk_graph, graph)) > 0) {
      order[psourceID] <- i
      wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% psourceID)
    }
    # All remaining outer nodes are to be tagged and pruned
    outerID <- names(which(degree(wk_graph, mode = "in") == 0))
    order[outerID] <- i
    wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% outerID)

    if (plot.it && vcount(wk_graph) > 0) {
      plot(wk_graph, vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1, main = i+1)
    }

    i <- i + 1
  }
  if (sum(is.na(order)) == 1) order[is.na(order)] <- i
  ### DONE

  order
}


#' @export
addOrder2Graph <- function(graph, order) {
  V(graph)$order <- ord
  V(graph)$color <- V(graph)$order

  # assign orders to edges
  for (i in 1:max(V(graph)$order)) {
    if (sum(V(graph)$order == i) > 1) {
      E(graph)[from(V(graph)[V(graph)$order == i])]$order <- i
    }
  }

  E(graph)$color <- rev(terrain.colors(max(E(graph)$order)))[E(graph)$order]
  E(graph)$width <- E(graph)$order / max(E(graph)$order) * 7

  graph
}



#' @export
addOrder2DRN <- function(rivs, graph) {
  rivs@data $ order <- E(graph)$order
  rivs@data $ width <- E(graph)$width
  rivs@data $ color <- E(graph)$color
  rivs@data $ start <- ends(g, E(g))[,1]
  rivs@data $ end <- ends(g, E(g))[,2]
  rivs <- addUpDownNodes(rivs)
  rivs
}


#' @export
findBuffer <- function(p, up_distance = 100, width = 25, searchWidth = NULL, demo = FALSE)
{
  message("finding xy shift... ")
  xyshift <- findShift(p, up_distance = 200)

  # calculate buffer using shifted sepa line
  # set sepa line (the search) buffer width based on river order
  message("finding buffer... ")

  # need to find search width first
  # find closest line segment on sepa network
  p_snap <- maptools::snapPointsToLines(p, rivs, 50)
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
  p_snap <- maptools::snapPointsToLines(p, rivs, 50)

  ## walk up sepa river
  out <- walkUpstream(p_snap, up_distance, useRiverOrder)
  p_upstr <- out$p_upstr
  seg3 <- out$seg

  # cut water polygon boundary here
  # tricky - one solution is
  seg3_shift <- shift(gLineMerge(seg3), shift[1], shift[2])
  buff_sepa <- rgeos::gBuffer(seg3_shift, width = sepaWidth, byid = FALSE, capStyle = "FLAT")

  cut_area <- rgeos::gIntersection(buff_sepa, wk_area, byid=TRUE, drop_lower_td=TRUE)
  cut_lines <- rgeos::gIntersection(buff_sepa, wk_lines, byid=TRUE, drop_lower_td=TRUE)

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
    
    out$cut_line <- gLineMerge(out$cut_line)
    out$cut_line@lines[[1]]@ID <- rownames(p@data)
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


#' @export
findShift <- function(p, up_distance = 100, useRiverOrder = TRUE) {
  # find the xy offset that minimises the distance between the sepa river line
  # and the river bank / line feature
  # need to first get a wide buffer idnoring river order
  # think about using down stream lines to help with alignment
  out <- getBuffer(p, up_distance = up_distance, width = 25, sepaWidth = 100, useRiverOrder = useRiverOrder)

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
makeDRNGraph <- function(rivs) {
  # Get edges from SEPA DRN and work out the nodes (top and tail the segments)
  updown_nodes <- getUpDownNodes(rivs)

  g <- make_graph(c(do.call(rbind, updown_nodes)), directed = TRUE)

  # add attributes to graph
  # location for nodes - this looks unnessisary, but its is safe.
  cc_vertices <- do.call(rbind, strsplit(V(g)$name, ":"))
  mode(cc_vertices) <- "numeric"
  V(g)$x <- cc_vertices[,1]
  V(g)$y <- cc_vertices[,2]

   # find mouths
  V(g)$color <- "blue"
  V(g)$color[degree(g, mode = "out") == 0] <- "red"

  # add lengths of edges
  edges <- lapply(rivs @ lines, function(x) x @ Lines[[1]] @ coords)
  E(g)$weight <- sapply(edges, LineLength, longlat = FALSE, sum = TRUE)
  g
}

#' @export
getUpDownNodes <- function(rivs) {
  edges <- lapply(rivs @ lines, function(x) x @ Lines[[1]] @ coords)
  cc_up_nodes <- t(sapply(edges, function(x) head(x,1)))
  cc_down_nodes <- t(sapply(edges, function(x) tail(x,1)))

  cc2label <- function(x) paste0(x, collapse = ":")
  up_nodes <- apply(cc_up_nodes, 1, cc2label)
  down_nodes <- apply(cc_down_nodes, 1, cc2label)

  list(up_node = up_nodes, down_node = down_nodes)
}

#' @export
addUpDownNodes <- function(rivs) {
  rivs@data <- cbind(rivs@data, getUpDownNodes(rivs))
  rivs
}


#' @export
summariseDRN <- function(g) {
  out <- list(
    mouths = V(g)[degree(g, mode = "out") == 0],
    sources = V(g)[degree(g, mode = "in") == 0],
    psuedonodes = V(g)[degree(g, mode = "in") == 1 & degree(g, mode = "out") == 1],
    confulences = V(g)[degree(g, mode = "in") == 2],
    braids = V(g)[degree(g, mode = "out") == 2],
    in3s = V(g)[degree(g, mode = "in") == 3],
    out3s = V(g)[degree(g, mode = "out") == 3],
    inout2plus = V(g)[degree(g, mode = "out") > 1 & degree(g, mode = "in") > 1],
    isdag = is.dag(g)
  )

  out
}



#' @export
openGoogleMaps <- function(point) {
  xy <- rowMeans(bbox(spTransform(point, CRS("+init=epsg:4326"))))
  url <- paste0("https://www.google.co.uk/maps/@",xy[2],",",xy[1],",3107m/data=!3m1!1e3")
  browseURL(url)
}


#' @export
plotbase <- function(xy, withGoogle = FALSE, plot.rivs = TRUE, lwd = 2, ...) {
  # crop things for quicker plotting and computation
  wk_area <- cropFeature(wareas, xy, buffer = 1000)
  wk_lines <- cropFeature(wlines, xy, buffer = 1000)
  if (plot.rivs) wk_rivs <- cropFeature(rivs, xy, buffer = 1000)

  if (withGoogle) {
    # transform to googlemap projection
    wk_area <- spTransform(wk_area, CRS(projection(r)))
    wk_lines <- spTransform(wk_lines, CRS(projection(r)))
    if (plot.rivs) wk_rivs <- spTransform(wk_rivs, CRS(projection(r)))

    r <- gmap(extent(spTransform(xy, CRS("+init=epsg:4326"))), type='satellite')
    plot(r, ...)
  }
  wk_area@bbox <- bbox(extent(xy))

  #
  plot(wk_area, add = withGoogle, border = NA, col = "lightblue", ...)
  if (plot.rivs) plot(wk_rivs, add = TRUE, col = "red", lwd = lwd)
  plot(wk_lines, add = TRUE, col = "blue", lwd = lwd)
}



#' @export
col_alpha <- function(col, alpha = 1) {
  alpha <- round(pmax(0, pmin(1, as.numeric(alpha))) * 255)
  alpha <- toupper(as.hexmode(alpha))
  paste0(sapply(col, function(x) colorRampPalette(x)(1)), alpha)
}
