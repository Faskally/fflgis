# Editted snap points to lines for river nodes

snapPointsToLinesbyOrder<-
  function (points, lines, maxDist = NA, withAttrs = TRUE, idField = NA)
  { # browser()
    if (rgeosStatus()) {
      if (!requireNamespace("rgeos", quietly = TRUE))
        stop("package rgeos required for snapPointsToLines")
    }
    else stop("rgeos not installed")
    if (class(points) == "SpatialPoints" && missing(withAttrs))
      withAttrs = FALSE
    if (class(points) == "SpatialPoints" && withAttrs == TRUE)
      stop("A SpatialPoints object has no attributes! Please set withAttrs as FALSE.")
    if (!is.na(maxDist)) {
      w = rgeos::gWithinDistance(points, lines, dist = maxDist,
                                 byid = TRUE)
      validPoints = apply(w, 2, any)
      validLines = apply(w, 1, any)
      points = points[validPoints, ]
      lines = lines[validLines, ]
    }
    # get a matrix of distances for the point to each of the lines (within the specified buffer area)
    d = rgeos::gDistance(points, lines, byid = TRUE)
    # Round to make sure you will always get a 0 because 0.00000000092766569 for example is
    # clearly 0 but won't come out in the which.
    d <- round(d, 2)
    # nearest_line_index = apply(d, 2, which.min) <<< ORIGINAL CODE
    # return all the rows in d which == 0 (rather than which.min that gives the first)
    # this gives you a vector of all the row positions of each river segment which has
    # 0 distance to the node
    all_nearest_line_indexes<-apply(d, 2, function(x) which(x==0))
    # now find out of all the closest line segments for each node with the highest order
    # have to wrap this in a a loop because which.max does not work if there is only 1 element
    # in the vector. Thus if there is only one line where there is 0 distance all_nearest_line_indexes
    # becomes nearest_line_index (although this should never happen apart from at the mouth as all other locations
    # will have at least 2 matches!)
    if(length(all_nearest_line_indexes) == 1) {
      nearest_line_index <- all_nearest_line_indexes
    } else {
      nearest_line_index<-apply(all_nearest_line_indexes, 2, function(x) x[which.max(lines[x, ]$order)])
    }
    # Get coordinates of all bits of the lines
    coordsLines = coordinates(lines)
    # Get coordinates of the points
    coordsPoints = coordinates(points)
    # Get coordinates of nearest line coordinate (for river nodes this will be the same as the node anyway)
    mNewCoords = vapply(1:length(points), function(x) nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                                                                         coordsPoints[x, ]), FUN.VALUE = c(0, 0))
    if (!is.na(idField))
      # ID of the line of the river line of interest in the lines sldf
      nearest_line_id = lines@data[, idField][nearest_line_index]
    else nearest_line_id = sapply(slot(lines, "lines"), function(i) slot(i,
                                                                         "ID"))[nearest_line_index]
    if (withAttrs)
      df = cbind(points@data, nearest_line_id)
    else df = data.frame(nearest_line_id, row.names = names(nearest_line_index))
    SpatialPointsDataFrame(coords = t(mNewCoords), data = df,
                           proj4string = CRS(proj4string(points)))
  }
