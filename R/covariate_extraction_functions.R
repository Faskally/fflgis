#' @export
extractGradientFromLine <- function(sl, dtm) {
  # line length
  length <- 
  # extract elevations
  cc <- coordinates(gLineMerge(sl))[[1]][[1]]
  p <- SpatialPoints(cc[c(1,nrow(cc)),], proj4string = crs(sl))
  el <- extract(dtm, p)
  # calculate gradient
  (el[1] - el[2]) / SpatialLinesLengths(sl)
}
