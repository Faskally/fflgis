
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
  int <- min(order(dists)[1:2])
  # new coords for line (cut upstream points out)
  cc_new <- rbind(coordinates(pt), cc[-(1:int),])
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

  rgeos::gIntersection(bbox, x, byid=TRUE, drop_lower_td=TRUE)
}



#' @export
getOrder <- function(graph, plot.it = FALSE) {

  stopifnot(is.named(graph))
  wk_graph <- graph
  
  # determine if the outermost nodes are to be pruned back
  idPsuedoSource <- function(gr, gr_full) {
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
  
  
  i <- 1
  order <- rep(NA, vcount(graph))
  names(order) <- V(graph)$name
  
  if (plot.it && vcount(wk_graph) > 0) {
    plot(wk_graph, vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1, main = i)
  }
  
  while(vcount(wk_graph) > 1) {
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
  #rivs@data $ start <- ends(g, E(g))[,1]
  #rivs@data $ end <- ends(g, E(g))[,2]
  rivs  
}




#' @export
getBuffer <- function(p, up_distance = 100, width = 25, sepaWidth = 50, shift = c(0,0), useRiverOrder = TRUE) {
  # expects to have the following objects available:
  #   rivs
  #   g
  #   wareas
  #   wlines
  #   dtm

  # have this here for now...
  require(maptools)
  require(deldir)

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
  #points(p_snap, col = "gold", pch = 16, cex = 2)

  ## move up sepa river
  # get segment that snapped point is on
  seg <- rivs[as.character(p_snap $ nearest_line_id),]
  # cut the segment to start at the snapped point
  seg2 <- cutLineDownstream(seg, p_snap)
  # how long is the line
  seg_len <- SpatialLinesLengths(seg2) # 17.4, so we will encouter a split
  if (seg_len < up_distance) {
    # what are the next segments?
    upVs <- names(V(g)[nei(seg@data[["start"]], mode = "in")])
    ids <- get.edge.ids(g, c(upVs[1], seg@data$start, upVs[2], seg@data$start))
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
  } else {
    # No need to cross junctions
    # get a point upstream
    p_upstr <- getPointUpstream(seg2, up_distance)
    points(p_upstr, col = "orange", pch = 16, cex = 2)
    # cut line upstream at new point
    seg3 <- cutLineUpstream(seg2, p_upstr)
  }

  #points(p_upstr, col = "orange", pch = 16, cex = 2)
  #plot(seg3, col = "brown", add = TRUE, lwd = 3)

  # cut water polygon boundary here
  # tricky - one solution is
  seg3_shift <- shift(gLineMerge(seg3), shift[1], shift[2])
  buff_sepa <- rgeos::gBuffer(seg3_shift, width = sepaWidth, byid = FALSE, capStyle = "FLAT")
  #plot(buff_sepa, col = paste0(grey(0.5), "22"), add = TRUE)

  cut_area <- rgeos::gIntersection(buff_sepa, wk_area, byid=TRUE, drop_lower_td=TRUE)
  cut_lines <- rgeos::gIntersection(buff_sepa, wk_lines, byid=TRUE, drop_lower_td=TRUE)

  # now add a buffer to the river segment
  buff_riva <- rgeos::gBuffer(cut_area, width = width, byid = FALSE)
  #plot(buff_riva, add = TRUE, col = paste0(grey(0.5), "77"), border = "green", lwd = 2)
  #plot(cut_area, col = "red", add = TRUE)

  # get buffer based on OS lines
  buff_rivl <- rgeos::gBuffer(cut_lines, width = width, byid = FALSE)
  #plot(buff_rivl, add = TRUE, col = paste0(grey(0.5), "77"), border = "green", lwd = 2)
  #plot(cut_lines, col = "pink", add = TRUE)

  # join the buffers incase there is a slight discrepancy
  buff_riv <- gUnion(buff_riva, buff_rivl)
  #plot(buff_riv, add = TRUE, col = paste0(grey(0.5), "77"), border = "cyan", lwd = 2)

  # return stuff
  list(buffer = buff_riv, p_upstr = p_upstr, riv_seg = seg3,
       cut_area = cut_area, cut_lines = cut_lines)  # gLineMerge(seg3)?
}


#' @export
findShift <- function(p, useRiverOrder = TRUE) {
  # find the xy offset that minimises the distance between the sepa river line
  # and the river bank / line feature
  # need to first get a wide buffer idnoring river order
  out <- getBuffer(p, up_distance = 100, width = 25, sepaWidth = 100, useRiverOrder = useRiverOrder)

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

  #  plotbase(xy, main = i)
  #  points(lxy)
  #  points(shift(lxy, par[1], par[2]))
  #  lines(shift(out$riv_seg, opt$par[1], opt$par[2]), lwd = 2)

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
