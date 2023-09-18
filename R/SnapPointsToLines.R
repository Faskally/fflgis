# New snapPointsToLines
#' @export

snapPointsToLines<-function (points, lines, maxDist = NA, withAttrs = TRUE, idField = NA)
{
  if (rgeosStatus()) {
    if (!requireNamespace("rgeos", quietly = TRUE))
      stop("package rgeos required for snapPointsToLines")
  }
  else stop("rgeos not installed")
  if (inherits(points, "SpatialPoints") && missing(withAttrs))
    withAttrs = FALSE
  if (inherits(points, "SpatialPoints") && withAttrs == TRUE)
    stop("A SpatialPoints object has no attributes! Please set withAttrs as FALSE.")
  d = gDistance(points, lines, byid = TRUE)
  if (!is.na(maxDist)) {
    distToLine <- apply(d, 2, min, na.rm = TRUE)
    validPoints <- distToLine <= maxDist
    distToPoint <- apply(d, 1, min, na.rm = TRUE)
    validLines <- distToPoint <= maxDist
    points <- points[validPoints, ]
    lines = lines[validLines, ]
    d = d[validLines, validPoints, drop = FALSE]
    distToLine <- distToLine[validPoints]
    if (!any(validPoints)) {
      if (is.na(idField)) {
        idCol = character(0)
      }
      else {
        idCol = lines@data[, idField][0]
      }
      newCols = data.frame(nearest_line_id = idCol, snap_dist = numeric(0))
      if (withAttrs)
        df <- cbind(points@data, newCols)
      else df <- newCols
      res <- SpatialPointsDataFrame(points, data = df,
                                    proj4string = CRS(proj4string(points)), match.ID = FALSE)
      return(res)
    }
  }
  else {
    distToLine = apply(d, 2, min, na.rm = TRUE)
  }
  nearest_line_index = apply(d, 2, which.min)
  coordsLines = coordinates(lines)
  coordsPoints = coordinates(points)
  mNewCoords = vapply(1:length(points), function(x) nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                                                                       coordsPoints[x, ]), FUN.VALUE = c(0, 0))
  if (!is.na(idField)) {
    nearest_line_id = lines@data[, idField][nearest_line_index]
  }
  else {
    nearest_line_id = sapply(slot(lines, "lines"), function(i) slot(i,
                                                                    "ID"))[nearest_line_index]
  }
  if (withAttrs)
    df = cbind(points@data, data.frame(nearest_line_id, snap_dist = distToLine))
  else df = data.frame(nearest_line_id, snap_dist = distToLine,
                       row.names = names(nearest_line_index))
  SpatialPointsDataFrame(coords = t(mNewCoords), data = df,
                         proj4string = CRS(proj4string(points)))
}
