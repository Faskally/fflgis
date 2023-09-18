#' A set of functions to help with
#'
#' @docType package
#' @name fflgis
#' @description Functions to help with spatial data analysis with a focus on river networks.
#' @importFrom igraph V E V<- E<- induced.subgraph igraph_opt igraph.options degree make_graph
#' @importFrom igraph vcount get.shortest.paths get.edge.ids is.dag is.named neighborhood
#' @importFrom sp 'spChFIDs<-'
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame SpatialLines
#' @importFrom sp SpatialLinesDataFrame SpatialLinesLengths
#' @importFrom sp SpatialPolygons SpatialPolygonsDataFrame Line LineLength Lines Polygon Polygons
#' @importFrom sp CRS spTransform proj4string spsample bbox is.projected
#' @importFrom raster projection extent crop shift extract crs 'crs<-' scalebar
#' @importFrom rgeos gBuffer gDistance gLineMerge gIntersection gUnion gUnaryUnion gIntersects gDifference gWithinDistance
#' @importFrom dismo gmap
#' @importFrom maptools snapPointsToLines rgeosStatus nearestPointOnLine
#' @importFrom rgdal readOGR writeOGR
#'
#' @importFrom grDevices colorRampPalette grey terrain.colors
#' @importFrom graphics lines points
#' @importFrom methods is slot
#' @importFrom stats optim
#' @importFrom utils browseURL combn head tail
NULL
