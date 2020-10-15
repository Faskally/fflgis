#' Find a point a given distance up stream on a river segment
#'
#' @param sldf a spatial line.
#' @param distance the distance to move upstream.
#' @return A spatial points object.
#' @export
getPointUpstream <- function(sldf, distance = 50) {
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
  else {
    df <- data.frame(Ind = Ind)
  }
  SpatialPointsDataFrame(midpoints, data = df, proj4string = crs(sldf))
}
