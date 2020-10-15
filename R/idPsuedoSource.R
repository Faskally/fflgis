# idPsuedoSource function

idPsuedoSource <- function(gr, gr_full) {



  # set local options for speed

  opt <- igraph_opt("return.vs.es")

  igraph.options(return.vs.es = FALSE)

  on.exit(igraph.options(return.vs.es = opt))



  # some usefull values

  deg.gr <- degree(gr, mode = "in")

  gr.names <- names(deg.gr)

  # determine if the outermost nodes are to be pruned back

  sourceID <- integer(0)



  # quick check for final main stem?

  if (sum(deg.gr == 0) == 1) {

    # then one mouth, one source and all nodes are to be removed

    gr.names[deg.gr == 0]

  }



  # otherwise get all psuedo nodes

  psuedoID <- which(deg.gr == 1)



  # identify if these nodes are one down from a source

  if (length(psuedoID) > 0) {

    # get upstream neighbours

    upID <- sapply(neighborhood(gr, 1, psuedoID, mode = "in"), function(x) x[2])

    # is the upID a source, if yes mark it to be removed

    sourceID <- c(sourceID, intersect(which(deg.gr == 0), upID))

  }



  # get all braided ends

  confluenceID <- which(deg.gr > 1)

  if (length(confluenceID) > 0) {

    # what are the upstream nodes

    #UPDATE R_igraph_neighborhood TO C_R_igraph_neighborhood..... 12/02/2020

    upID <- lapply(.Call(igraph:::C_R_igraph_neighborhood, gr,

                         igraph:::as.igraph.vs(gr, confluenceID) - 1, as.numeric(1),

                         as.numeric(2), as.integer(0), PACKAGE = "igraph"),

                   function(x) x[-1] + 1)



    # are these upstream nodes of confluences sources?

    is.PsuedoSource <- sapply(upID, function(id) sum(deg.gr[id]) == 0)

    # if so are these upstream nodes connected?

    if (any(is.PsuedoSource)) {

      upIDps <- upID[is.PsuedoSource]

      len <- sapply(upIDps, length)

      is.braid <- rep(FALSE, length(len))

      # odd one... why would this happen??

      if (any(len == 1)) {

        is.braid[len == 1] <- TRUE

      }



      gr.names.map <- as.vector(V(gr_full)[gr.names])

      if (any(len == 2)) {

        is.braid[len == 2] <- sapply(upIDps[len == 2], function(x) {

          upcommon <- do.call(intersect,

                              .Call(igraph:::C_R_igraph_neighborhood, gr_full,

                                    gr.names.map[x] - 1, as.numeric(20),

                                    as.numeric(2), as.integer(0), PACKAGE = "igraph"))

          length(upcommon) > 0

        })

      }



      if (any(len > 2)) {

        is.braid[len > 2] <- sapply(upIDps[len > 2], function(x) {

          upIDs <- .Call(igraph:::C_R_igraph_neighborhood, gr_full,

                         gr.names.map[x] - 1, as.numeric(20),

                         as.numeric(2), as.integer(0), PACKAGE = "igraph")

          # if there is atleast one unconnected node then this is

          # a confulence,

          # otherwise if all up nodes are connected then it is a braid.

          upcommon <- apply(combn(1:length(upIDs[c(1,1,2)]), 2), 2, function(ids) length(do.call(intersect, upIDs[ids])))

          all(upcommon > 0)

        })

      }

      # mark confulenceIDs that are braids as sources to be removed

      sourceID <- c(sourceID, unlist(upIDps[is.braid]))

    }

  }



  gr.names[sourceID]

}
