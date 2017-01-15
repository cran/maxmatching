#' Transforming a matrix representing a bipartite graph to an igraph
#' @title bipartite_matrix_to_igraph
#' @name bipartite_matrix_to_igraph
#' @description Transform a bipartite matrix to an igraph
#' @details Transform a bipartite matrix to an igraph
#' @param bipartite_matrix A matrix with rows representing group1,
#' and columns representing group 2,
#' the element at row i column j is the weight of the edge from i-th vertex in group 1
#' to j-th vertex in group 2
#' @param weighted True if the bipartite_matrix has weight, false otherwise
#'
#' @return The igraph
bipartite_matrix_to_igraph <- function(bipartite_matrix, weighted = FALSE) {
  n1 <- nrow(bipartite_matrix)
  n2 <- ncol(bipartite_matrix)
  group1_id <- 1:n1
  group2_id <- n1 + 1:n2
  row_id <- c(outer(group1_id, group2_id, FUN = function(x, y) x))
  col_id <- c(outer(group1_id, group2_id, FUN = function(x, y) y))
  weights <- c(bipartite_matrix)
  if (weighted) {
    edges <- c(rbind(row_id, col_id))
    G <- igraph::graph(edges)
    igraph::E(G)$weight <- weights
  } else {
    row_id <- row_id[weights == 1]
    col_id <- col_id[weights == 1]
    edges <- c(rbind(row_id, col_id))
    G <- igraph::graph(edges)
  }
  return(G)
}

# Interface to test the blossom algorithm
blossom_test <- function(edges, weights, maxcardinality=FALSE) {
	G <- igraph::graph(edges, directed=FALSE)
	igraph::E(G)$weight <- weights
	return(blossom(G, weighted=TRUE, maxcardinality=maxcardinality)$mate)
}


testthat::test_that("hungarian_test1", {
  x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3)
  hungarian_G1 <- bipartite_matrix_to_igraph(x, weighted = TRUE)
  matching_res <- maxmatching(hungarian_G1, weighted = TRUE)
  testthat::expect_equal(matching_res$matching_weight, 14)
})

testthat::test_that("hungarian_test2", {
  n1 <- 20
  n2 <- 20
  x <- matrix(runif(n1 * n2), nrow = n1, ncol = n2)
  hungarian_G2 <- bipartite_matrix_to_igraph(x, weighted = TRUE)
  matching_res <- maxmatching(hungarian_G2, weighted = TRUE)
  require(clue)
  lsap_assignment <- clue::solve_LSAP(x, maximum = TRUE)
  lsap_weight <- sum(x[cbind(seq_along(lsap_assignment), lsap_assignment)])
  testthat::expect_equal(matching_res$matching_weight, lsap_weight)
})

testthat::test_that("bipartite_matching_test1", {
  n1 <- 20
  n2 <- 20
  x <- matrix(sample(c(0, 1), n1 * n2, replace = TRUE), nrow = n1, ncol = n2)
  bipartite_G1 <- bipartite_matrix_to_igraph(x, weighted = FALSE)
  matching_res <- maxmatching(bipartite_G1, weighted = FALSE)
  lsap_assignment <- clue::solve_LSAP(x, maximum = TRUE)
  lsap_size <- sum(x[cbind(seq_along(lsap_assignment), lsap_assignment)])
  testthat::expect_equal(matching_res$matching_size, lsap_size)
})

# TODO: Add test cases here
testthat::test_that("blossom_test", {
  # single edge
  expect_equal(blossom_test(c(1,2), c(1)), c(2,1))
  expect_equal(blossom_test(c(1,2,2,3), c(10,11)), c(-1,3,2))
  expect_equal(blossom_test(c(1,2,2,3,3,4), c(5,11,5)), c(-1,3,2,-1))
  # maximum cardinality
  expect_equal(blossom_test(c(1,2,2,3,3,4), c(5,11,5), TRUE), c(2,1,4,3))
  # floating point weights
  expect_equal(blossom_test(c(1,2,2,3,1,3,1,4), c(pi,exp(1),3.0,sqrt(2.0))), c(4,3,2,1))
  # negative weights
  expect_equal(blossom_test(c(1,2,1,3,2,3,2,4,3,4), c(2,-2,1,-1,-6)), c(2,1,-1,-1))
  expect_equal(blossom_test(c(1,2,1,3,2,3,2,4,3,4), c(2,-2,1,-1,-6), TRUE), c(3,4,1,2))
  # create S-blossom and use it for augmentation
  expect_equal(blossom_test(c(1,2,1,3,2,3,3,4), c(8,9,10,7)), c(2,1,4,3))
  expect_equal(blossom_test(c(1,2,1,3,2,3,3,4,1,6,4,5), c(8,9,10,7,5,6)), c(6,3,2,5,4,1))
  # create S-blossom, relabel as T-blossom, use for augmentation
  expect_equal(blossom_test(c(1,2,1,3,2,3,1,4,4,5,1,6), c(9,8,10,5,4,3)), c(6,3,2,5,4,1))
  expect_equal(blossom_test(c(1,2,1,3,2,3,1,4,4,5,1,6), c(9,8,10,5,3,4)), c(6,3,2,5,4,1))
  expect_equal(blossom_test(c(1,2,1,3,2,3,1,4,4,5,3,6), c(9,8,10,5,3,4)), c(2,1,6,5,4,3))
  # create nested S-blossom, use for augmentation
  expect_equal(blossom_test(c(1,2,1,3,2,3,2,4,3,5,4,5,5,6), c(9,9,10,8,8,10,6)), c(3,4,1,2,6,5))
  # create S-blossom, relabel as S, include in nested S-blossom
  expect_equal(blossom_test(c(1,2,1,7,2,3,3,4,3,5,4,5,5,6,6,7,7,8), c(10,10,12,20,20,25,10,10,8)), c(2,1,4,3,6,5,8,7))
  # create nested S-blossom, augment, expand recursively
  expect_equal(blossom_test(c(1,2,1,3,2,3,2,4,3,5,4,5,4,6,5,7,6,7,7,8), c(8,8,10,12,12,14,12,12,14,12)), c(2,1,5,6,3,4,8,7))
  # create S-blossom, relabel as T, expand
  expect_equal(blossom_test(c(1,2,1,5,1,6,2,3,3,4,4,5,4,8,5,7), c(23,22,15,25,22,25,14,13)), c(6,3,2,8,7,1,5,4))
  # create nested S-blossom, relabel as T, expand
  expect_equal(blossom_test(c(1,2,1,3,1,8,2,3,2,4,3,5,4,5,4,7,5,6), c(19,20,8,25,18,18,13,7,7)), c(8,3,2,7,6,5,4,1))
  # create blossom, relabel as T in more than one way, expand, augment
  expect_equal(blossom_test(c(1,2,1,5,2,3,3,4,4,5,1,6,3,9,4,8,5,7,9,10), c(45,45,50,45,50,30,35,35,26,5)), c(6,3,2,8,7,1,5,4,10,9))
  # again but slightly different
  expect_equal(blossom_test(c(1,2,1,5,2,3,3,4,4,5,1,6,3,9,4,8,5,7,9,10), c(45,45,50,45,50,30,35,26,40,5)), c(6,3,2,8,7,1,5,4,10,9))
  # create blossom, relabel as T, expand such that a new least-slack S-to-free edge is produced, augment
  expect_equal(blossom_test(c(1,2,1,5,2,3,3,4,4,5,1,6,3,9,4,8,5,7,9,10), c(45,45,50,45,50,30,35,28,26,5)), c(6,3,2,8,7,1,5,4,10,9))
  # create nested blossom, relabel as T in more than one way, expand outer blossom such that inner blossom ends up on an augmenting path
  expect_equal(blossom_test(c(1,2,1,7,2,3,3,4,4,5,4,6,5,6,6,7,1,8,3,11,5,9,7,10,11,12), c(45,45,50,45,95,94,94,50,30,35,36,26,5)), c(8,3,2,6,9,4,10,1,5,7,12,11))
  # create nested S-blossom, relabel as S, expand recursively
  expect_equal(blossom_test(c(1,2,1,3,2,3,2,4,3,5,4,5,1,8,5,7,7,6,8,10,4,9), c(40,40,60,55,55,50,15,30,10,10,30)), c(2,1,5,9,3,7,6,10,4,8))
})
