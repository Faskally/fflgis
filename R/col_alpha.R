#' @export
col_alpha <- function(col, alpha = 1) {
  alpha <- round(pmax(0, pmin(1, as.numeric(alpha))) * 255)
  alpha <- toupper(as.hexmode(alpha))
  paste0(sapply(col, function(x) colorRampPalette(x)(1)), alpha)
}
