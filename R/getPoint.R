#' @export
getPoint <- function(coords, distance) {
  coords <- coords[nrow(coords):1, ]
  dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
  dist_mid <- distance
  if (sum(dist) < distance) {
    return(coords[nrow(coords), ])
  }
  dist_cum <- c(0, cumsum(dist))
  end_index <- which(dist_cum > dist_mid)[1]
  start_index <- end_index - 1
  start <- coords[start_index, ]
  end <- coords[end_index, ]
  dist_remaining <- dist_mid - dist_cum[start_index]
  mid <- start + (end - start) * (dist_remaining / dist[start_index])
  return(mid)
}
