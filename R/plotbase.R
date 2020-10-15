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

    r <- gmap(extent(spTransform(xy, CRS("+init=epsg:4326"))), type = "satellite")
    plot(r, ...)
  }
  wk_area@bbox <- bbox(extent(xy))

  #
  plot(wk_area, add = withGoogle, border = NA, col = "lightblue", ...)
  if (plot.rivs) plot(wk_rivs, add = TRUE, col = "red", lwd = lwd)
  plot(wk_lines, add = TRUE, col = "blue", lwd = lwd)
}
