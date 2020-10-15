#' @export
openGoogleMaps <- function(point) {
  xy <- rowMeans(bbox(spTransform(point, CRS("+init=epsg:4326"))))
  url <- paste0("https://www.google.co.uk/maps/@", xy[2], ",", xy[1], ",3107m/data=!3m1!1e3")
  browseURL(url)
}
