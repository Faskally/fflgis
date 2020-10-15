#' Cut off the upstream segment of a river segment at a point
#'
#' @param line a spatial line.
#' @param pt a spatial point, previously snapped to \code{line}.
#' @return A spatial lines object.
#' @examples
#' bng <- CRS("+init=epsg:27700")
#' cc <- cbind(c(0,0), c(0,1))
#' line <- SpatialLines(list(Lines(list(Line(cc)), ID = "A")), bng)
#' pt <- SpatialPoints(cbind(0, .5), bng)
#' up_line <- cutLineUpstream(line, pt)
#'
#' plot(line)
#' points(cc)
#' points(pt, col = "green")
#' lines(up_line, col = "red")
#' y <- function(x) 0.1*x^2
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
cutLineUpstream <- function(line, pt) {
  # get digitised points along river
  cc <- coordinates(line)[[1]][[1]]
  # find the distance between the snapped point and the digitised points
  dists <- sqrt(colSums((t(cc) - c(coordinates(pt)))^2))
  # check if the point lies exactly on the first point
  if (dists[length(dists)] == 0) {
    warning("line is completely cut, i.e. pt is same as last point in line")
  }
  # find which interval the snapped point is in
  # split into many line segments
  slcc <- SpatialLines(lapply(1:(nrow(cc)-1), function(i) Lines(list(Line(cc[i + 0:1,])), ID = i)), crs(line))
  d <- rgeos::gDistance(pt, slcc, byid = TRUE)
  int <- which.min(d)
  # new coords for line (cut upstream points out)
  cc_new <- rbind(coordinates(pt), cc[(int+1):nrow(cc),])
  SpatialLines(list(Lines(list(Line(cc_new)), ID = "A")), crs(line))
}
