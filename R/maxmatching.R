#' maxmatching
#' @title Maximum Matching
#' @name maxmatching
#' @description Compute the maximum matching for undirected graph
#' @details TODO
#' @param G undirected igraph object representing the input
#' @param weighted whether the graph is weighted, if the graph is weighted, the weight should be able to be accessed with igraph::E(G)$weight
#' @param maxcardinality Ignore if the graph is bipartite, and unmeaningful if the graph is unweighted. If the graph is non-bipartite and weighted, only computes the maximum weighted matching among all maximum cardinality matchings.
#'
#' @return The matchings in a list
#'
#' @examples
#' # Unweighted general graph
#' G1 <- igraph::graph(c(1, 2, 1, 3, 1, 4, 3, 4, 3, 5,
#' 5, 6, 6, 7, 7, 8, 8, 9, 3, 8, 5, 8), directed = FALSE)
#' maxmatching(G1, weighted = FALSE)
#' # Unweighted bipartite graph
#' G2 <- igraph::graph(c(1, 5, 1, 6, 1, 7, 2, 5, 2, 8,
#' 3, 6, 3, 7, 3, 8, 4, 6, 4, 7, 4, 8), directed = FALSE)
#' maxmatching(G2, weighted = FALSE)
#' # Weighted bipartite graph
#' G3 <- igraph::graph(c(1, 5, 1, 6, 1, 7, 2, 5, 2, 8,
#' 3, 6, 3, 7, 3, 8, 4, 6, 4, 7, 4, 8), directed = FALSE)
#' igraph::E(G3)$weight <- runif(length(igraph::E(G3)), 0, 1)
#' maxmatching(G3, weighted = TRUE)
#'
#'
#' @import igraph
#' @export
maxmatching <- function(G, weighted=FALSE,maxcardinality=FALSE) {
  bipartite_matching <- igraph::bipartite.mapping(G)
  if (bipartite_matching$res) {
    # Source code for igraph matching
    # https://github.com/igraph/igraph/blob/02005680aebb4281f2964c6d2d8d98851fce63b9/src/matching.c
    if (weighted) {
      return(igraph::max_bipartite_match(graph = G, types = bipartite_matching$type, weights = igraph::E(G)$weight))
    } else {
      return(igraph::max_bipartite_match(graph = G, types = bipartite_matching$type))
    }
  } else {
    return(blossom(G, weighted=weighted, maxcardinality=maxcardinality))
  }
}
