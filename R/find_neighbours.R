#' Find matrix of n_neigh neighbours for each node
#'
#' This function computes the matrix of n_neigh nearest neighbors for each node based on their coordinates (nodes_x and nodes_y).
#' The result is a matrix where each row corresponds to a node, and the columns represent its n_neigh nearest neighbors.
#' 
#' @param nodes_x A numeric vector containing the x-coordinates of the nodes.
#' @param nodes_y A numeric vector containing the y-coordinates of the nodes.
#' @param n_neigh An integer specifying the number of nearest neighbors to find for each node (excluding itself).
#' @return A matrix where each row corresponds to a node, and the columns represent its n_neigh nearest neighbors.
#' @import FNN
#' @export
#' 
#' @examples
#' nodes_x <- c(1, 2, 3, 4, 5)
#' nodes_y <- c(2, 4, 6, 8, 10)
#' n_neigh <- 2
#' find_neighbours(nodes_x, nodes_y, n_neigh)
find_neighbours <- function(nodes_x, nodes_y, n_neigh) {
  
  library(FNN)
  
  # Function implementation
  n <- length(nodes_x)
  
  nodes <- data.frame(x = nodes_x, y = nodes_y)
  
  # Find n_neigh nearest neighbors
  lst_neighbours <- FNN::get.knn(data = nodes, k = n_neigh - 1)$nn.index
  
  lst_neighbours <- cbind(1:n, lst_neighbours)
  
  return(lst_neighbours)
}