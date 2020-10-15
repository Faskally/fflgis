#' @export
summariseDRN <- function(g) {
  out <- list(
    mouths = V(g)[degree(g, mode = "out") == 0],
    sources = V(g)[degree(g, mode = "in") == 0],
    psuedonodes = V(g)[degree(g, mode = "in") == 1 & degree(g, mode = "out") == 1],
    confulences = V(g)[degree(g, mode = "in") == 2],
    braids = V(g)[degree(g, mode = "out") == 2],
    in3s = V(g)[degree(g, mode = "in") == 3],
    out3s = V(g)[degree(g, mode = "out") == 3],
    inout2plus = V(g)[degree(g, mode = "out") > 1 & degree(g, mode = "in") > 1],
    isdag = is.dag(g)
  )

  out
}
