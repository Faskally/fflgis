#' Cut off the down stream segment of a river segment at a point
#'
#' @param line a spatial line.
#' @param pt a spatial point, previously snapped to \code{line}.
#' @return A spatial lines object.
#' @examples
#' library(sp)
#' bng <- CRS("+init=epsg:27700")
#' cc <- cbind(c(0,0), c(0,1))
#' line <- SpatialLines(list(Lines(list(Line(cc)), ID = "A")), bng)
#' pt <- SpatialPoints(cbind(0, .5), bng)
#' down_line <- cutLineDownstream(line, pt)
#'
#' plot(line)
#' points(cc)
#' points(pt, col = "green")
#' lines(down_line, col = "red")
#' y <- function(x) 0.1*x^2
#'
#' cc <- cbind(-10:10, y(-10:10))
#' line <- SpatialLines(list(Lines(list(Line(cc)), ID = "A")), bng)
#' pt <- SpatialPoints(cbind(-8.020474, 6.493534), bng)
#' pt <- snapPointsToLines(pt, line)
#' down_line <- cutLineDownstream(line, pt)
#'
#' plot(down_line, col = "red", lwd = 2, xlim = c(-10, -2), ylim = c(5, 10))
#' lines(line)
#' points(cc)
#' points(pt, col = "green")
#' @export
cutLineDownstream <- function(line, pt) {
  # get digitised points along river
  cc <- coordinates(line)[[1]][[1]]
  # find the distances between the snapped point
  # and the digitised points
  dists <- sqrt(colSums((t(cc) - c(coordinates(pt)))^2))
  # check if the point lies exactly on the first point
  if (dists[1] == 0) {
    warning("line is completely cut, i.e. pt is same as first point in line")
  }
  # find which interval the snapped point is in
  # split into many line segments
  slcc <- SpatialLines(lapply(1:(nrow(cc)-1), function(i) Lines(list(Line(cc[i + 0:1,])), ID = i)), crs(line))
  d <- gDistance(pt, slcc, byid = TRUE)
  int <- which.min(d)
  # new coords for line (cut downstream points out)
  cc_new <- rbind(cc[1:int,],coordinates(pt))
  SpatialLines(list(Lines(list(Line(cc_new)), ID = "A")), crs(line))
}
