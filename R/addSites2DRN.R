#' @export
addSites2DRN <- function(sites, rivs, site_names, mdist = 200) {
  message("Make sure that sites locations are unique if possible. \nOtherwise site names could be overwritten.\n")
  if (is.null(names(site_names))) {
    rname <- "site.name"
  }
  else {
    rname <- names(site_names)
  }
  if (is.null(rivs[[rname]])) {
    rivs[[rname]] <- ""
  }
  site_names <- sites[[site_names]]
  for (i in 1:length(sites)) {
    rivs <- addSite2DRN(sites[i, ], rivs, site_names[i],
      rname,
      mdist = mdist
    )
  }
  rivs
}
