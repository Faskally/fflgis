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
  cc_new <- rbind(cc[1:(int-1),],coordinates(pt))
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


#' Find a point a given distance up stream on a river segment
#'
#' @param line a spatial line.
#' @param dist the distance to move upstream.
#' @return A spatial lines object.
#' @examples
#' bng <- CRS("+init=epsg:27700")
#' cc <- cbind(c(0,0), c(0,1))
#' line <- SpatialLines(list(Lines(list(Line(cc)), ID = "A")), bng)
#' up_pt <- getPointUpstream(line, 0.5)
#'
#' plot(line)
#' points(cc)
#' points(up_pt, col = "green")
#'
#' cc <- cbind(-10:10, y(-10:10))
#' line <- SpatialLines(list(Lines(list(Line(cc)), ID = "A")), bng)
#' pt <- SpatialPoints(cbind(-8.020474, 6.493534), bng)
#' pt <- snapPointsToLines(pt, line)
#' up_line <- cutLineUpstream(line, pt)
#'
#' plot(up_line, col = "red", lwd = 2, xlim = c(-10, -2), ylim = c(5, 10))
#' lines(line)
#' points(cc)
#' points(pt, col = "green")
#' @export
getPointUpstream <- function (sldf, distance = 50)
{
  stopifnot(is.projected(sldf))
  Lns <- slot(sldf, "lines")
  hash_lns <- sapply(Lns, function(x) length(slot(x, "Lines")))
  N <- sum(hash_lns)
  midpoints <- matrix(NA, ncol = 2, nrow = N)
  Ind <- integer(length = N)
  ii <- 1
  for (i in 1:length(Lns)) {
    Lnsi <- slot(Lns[[i]], "Lines")
    for (j in 1:hash_lns[i]) {
      Ind[ii] <- i
      midpoints[ii, ] <- getPoint(slot(Lnsi[[j]], "coords"), distance = distance)
      ii <- ii + 1
    }
  }
  if (is(sldf, "SpatialLinesDataFrame")) {
    df0 <- slot(sldf, "data")[Ind, ]
    df <- as.data.frame(cbind(df0, Ind))
  }
  else df <- data.frame(Ind = Ind)
  SpatialPointsDataFrame(midpoints, data = df, proj4string = crs(sldf))
}


#' @export
getPoint <- function(coords, distance)
{
  coords <- coords[nrow(coords):1,]
  dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
  dist_mid <- distance
  if (sum(dist) < distance) {
    return(coords[nrow(coords),])
  }
  dist_cum <- c(0, cumsum(dist))
  end_index <- which(dist_cum > dist_mid)[1]
  start_index <- end_index - 1
  start <- coords[start_index, ]
  end <- coords[end_index, ]
  dist_remaining <- dist_mid - dist_cum[start_index]
  mid <- start + (end - start) * (dist_remaining/dist[start_index])
  return(mid)
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
  upDownNodes <- getUpDownNodes(rivs)
  rivs$start <- upDownNodes$up_node
  rivs$end <- upDownNodes$down_node
  rivs$up_node <- upDownNodes$up_node
  rivs$down_node <- upDownNodes$down_node
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
plotbase <- function(xy, withGoogle = FALSE, plot.rivs = TRUE, lwd = 2, wareas, wlines, rivs, ...) {
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



# ---------------------------------------------------------------------------------------------------
# added from seaDistance
# ---------------------------------------------------------------------------------------------------

#' @export
getSeaDistance <- function(p, rivs, g) {
  p_snap <- maptools::snapPointsToLines(p, rivs)

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



#' @export
plotPath <- function(path, g) {
  vgraph <- induced.subgraph(g, path)
  E(vgraph)$width <- 3
  E(vgraph)$color <- "orange"
  E(vgraph)$arrow.size <- 0.1
  V(vgraph)$label <- NA
  V(vgraph)$size <- 0.1

  plot(vgraph, add=TRUE, rescale=FALSE)
}



#' @export
idPsuedoSource <- function(gr, gr_full) {

  # set local options for speed
  opt <- igraph_opt("return.vs.es")
  igraph.options(return.vs.es = FALSE)
  on.exit(igraph.options(return.vs.es = opt))

  # some usefull values
  deg.gr <- degree(gr, mode = "in")
  gr.names <- names(deg.gr)
  # determine if the outermost nodes are to be pruned back
  sourceID <- integer(0)

  # quick check for final main stem?
  if (sum(deg.gr == 0) == 1) {
    # then one mouth, one source and all nodes are to be removed
    gr.names[deg.gr == 0]
  }

  # otherwise get all psuedo nodes
  psuedoID <- which(deg.gr == 1)

  # identify if these nodes are one down from a source
  if (length(psuedoID) > 0) {
    # get upstream neighbours
    upID <- sapply(neighborhood(gr, 1, psuedoID, mode = "in"), function(x) x[2])
    # is the upID a source, if yes mark it to be removed
    sourceID <- c(sourceID, intersect(which(deg.gr == 0), upID))
  }

  # get all braided ends
  confluenceID <- which(deg.gr > 1)
  if (length(confluenceID) > 0) {
    # what are the upstream nodes
    upID <- lapply(.Call("R_igraph_neighborhood", gr,
                 igraph:::as.igraph.vs(gr, confluenceID) - 1, as.numeric(1),
                 as.numeric(2), as.integer(0), PACKAGE = "igraph"),
            function(x) x[-1] + 1)

    # are these upstream nodes of confluences sources?
    is.PsuedoSource <- sapply(upID, function(id) sum(deg.gr[id]) == 0)
    # if so are these upstream nodes connected?
    if (any(is.PsuedoSource)) {
      upIDps <- upID[is.PsuedoSource]
      len <- sapply(upIDps, length)
      is.braid <- rep(FALSE, length(len))
      # odd one... why would this happen??
      if (any(len == 1)) {
        is.braid[len == 1] <- TRUE
      }

      gr.names.map <- as.vector(V(gr_full)[gr.names])
      if (any(len == 2)) {
        is.braid[len == 2] <- sapply(upIDps[len == 2], function(x) {
          upcommon <- do.call(intersect,
                              .Call("R_igraph_neighborhood", gr_full,
                                    gr.names.map[x] - 1, as.numeric(20),
                                    as.numeric(2), as.integer(0), PACKAGE = "igraph"))
          length(upcommon) > 0
        })
      }

      if (any(len > 2)) {
        is.braid[len > 2] <- sapply(upIDps[len > 2], function(x) {
          upIDs <- .Call("R_igraph_neighborhood", gr_full,
                         gr.names.map[x] - 1, as.numeric(20),
                         as.numeric(2), as.integer(0), PACKAGE = "igraph")
          # if there is atleast one unconnected node then this is
          # a confulence,
          # otherwise if all up nodes are connected then it is a braid.
          upcommon <- apply(combn(1:length(upIDs[c(1,1,2)]), 2), 2, function(ids) length(do.call(intersect, upIDs[ids])))
          all(upcommon > 0)
          })
      }
      # mark confulenceIDs that are braids as sources to be removed
      sourceID <- c(sourceID, unlist(upIDps[is.braid]))
    }
  }

  gr.names[sourceID]
}


#' @export
stripOrder <- function(graph) {

  # set local options
  opt <- igraph_opt("return.vs.es")
  igraph.options(return.vs.es = TRUE)
  on.exit(igraph.options(return.vs.es = opt))

  stopifnot(is.named(graph))
  wk_graph <- graph

  # lop while there are psuedo source nodes next to sources
  while(length(psourceID <- idPsuedoSource(wk_graph, graph)) > 0) {
    #print(psourceID)
    wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% psourceID)
  }
  # All remaining outer nodes are to be tagged and pruned
  outerID <- names(which(degree(wk_graph, mode = "in") == 0))
  if (length(outerID)==0) {
    message("Possible problem with graph: no source nodes but nvertices > 1\n",
            "This can happen if there is a cycle in the graph, or an edge pointing backwards.")
    return(wk_graph)
  }
  wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% outerID)

  wk_graph
}



#' @export
getOrder <- function(graph, plot.it = FALSE, max_order = 10) {

  # set local options
  opt <- igraph_opt("return.vs.es")
  igraph.options(return.vs.es = TRUE)
  on.exit(igraph.options(return.vs.es = opt))

  stopifnot(is.named(graph))
  wk_graph <- graph

  i <- 1
  order <- rep(NA, vcount(graph))
  names(order) <- V(graph)$name

  if (plot.it && vcount(wk_graph) > 0) {
    plot(wk_graph, vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1,
         main = i, rescale = FALSE,
         xlim = range(V(graph)$x), ylim = range(V(graph)$y))
  }

  while(vcount(wk_graph) > 1) {
    if (i > max_order) {
      message("Possible problem with graph: maximum order reached.  \n",
              "Returning graph.")
      return(wk_graph)
    }

    #print(i)
    # lop while there are psuedo source nodes next to sources
    while(length(psourceID <- idPsuedoSource(wk_graph, graph)) > 0) {
      #print(psourceID)
      order[psourceID] <- i
      wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% psourceID)
    }
    # All remaining outer nodes are to be tagged and pruned
    outerID <- names(which(degree(wk_graph, mode = "in") == 0))
    if (length(outerID)==0) {
      message("Possible problem with graph: no source nodes but nvertices > 1\n",
               "This can happen if there is a cycle in the graph, or an edge pointing backwards.\n",
               "Returning graph.")
      return(wk_graph)
    } #else if (length(outerID)==1) {
      # then we are on the max order?
    #}
    order[outerID] <- i
    wk_graph <- induced.subgraph(wk_graph, !V(wk_graph)$name %in% outerID)

    if (plot.it && vcount(wk_graph) > 0) {
      plot(wk_graph, vertex.label = NA, vertex.size = 0, edge.arrow.size = 0.1,
           main = i+1, rescale = FALSE,
           xlim = range(V(graph)$x), ylim = range(V(graph)$y))
    }

    i <- i + 1
  }
  if (sum(is.na(order)) == 1) order[is.na(order)] <- i
  ### DONE

  order
}


#' @export
addOrder2Graph <- function(graph, order) {
  V(graph)$order <- order
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
  #rivs@data $ start <- ends(graph, E(graph))[,1]
  #rivs@data $ end <- ends(graph, E(graph))[,2]
  rivs <- addUpDownNodes(rivs)
  rivs
}




addSite2DRN <- function(site, rivs, site_name, rname) {
  # find segment to split
  snap_site <- snapPointsToLines(site, rivs)

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



#' @export
addSites2DRN <- function(sites, rivs, site_names) {
  message("Make sure that sites locations are unique if possible. \nOtherwise site names could be overwritten.\n")
  # get column name to put site names in
  if (is.null(names(site_names))) {
    rname <- "site.name"
  } else {
    rname <- names(site_names)
  }
  # if rivs does not have site columns in dataframe, then add them as ""
  if (is.null(rivs[[rname]])) {
    rivs[[rname]] <- ""
  }
  # get site name from sites data
  site_names <- sites[[site_names]]
  for (i in 1:length(sites)) {
    rivs <- addSite2DRN(sites[i,], rivs, site_names[i], rname)
  }
  rivs
}




