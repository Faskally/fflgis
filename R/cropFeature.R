#' @export
cropFeature <- function(x, xy, buffer = 0) {
  # make a spatial polygon for the new limits
  bbox <- cbind(
    x = sort(xy$x)[c(1, 2, 2, 1, 1)],
    y = sort(xy$y)[c(1, 1, 2, 2, 1)]
  )
  bbox <- SpatialPolygons(list(Polygons(list(Polygon(bbox)), ID = "a")))
  bbox @ proj4string <- x @ proj4string

  if (buffer > 0) {
    bbox <- gBuffer(bbox, width = buffer, byid = TRUE)
  }

  ind <- as.vector(rgeos::gIntersects(x, bbox, byid = TRUE))

  if (sum(ind) == 0) ind <- 1

  x[ind, ]
}
