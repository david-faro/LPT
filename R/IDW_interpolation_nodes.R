#' Inverse Distance Weighted (IDW) Interpolation for Flow Field Data
#'
#' This function performs Inverse Distance Weighted (IDW) interpolation to estimate the flow velocities (x and y components) and turbulent kinematic properties (K) at a given point (x, y) in a flow field based on the values at neighboring nodes. IDW is a simple spatial interpolation method that gives higher weight to closer neighboring nodes, resulting in a smooth interpolated value.
#'
#' @param x The x-coordinate of the point where the interpolation is performed.
#' @param y The y-coordinate of the point where the interpolation is performed.
#' @param nodes_x A numeric vector containing the x-coordinates of the neighboring nodes.
#' @param nodes_y A numeric vector containing the y-coordinates of the neighboring nodes.
#' @param vel_x A numeric vector containing the x-component of flow velocities at the neighboring nodes.
#' @param vel_y A numeric vector containing the y-component of flow velocities at the neighboring nodes.
#' @param K A numeric vector containing the turbulent kinematic properties (K) at the neighboring nodes.
#' @param id.node The index of the node for which the interpolation is being performed.
#' @param lst_neighbours A matrix representing the matrix of nearest neighbors for each node in the flow field.
#' @return A data frame containing the interpolated flow velocities (vel_x and vel_y) and the turbulent kinematic property (int_K) at the given point (x, y). The data frame also includes the index of the closest node (node.id) used for the interpolation.
#' @export
#'
#' @examples
#' # Sample flow field data
#' nodes_x <- c(1, 2, 3, 4, 5)
#' nodes_y <- c(2, 4, 6, 8, 10)
#' vel_x <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' vel_y <- c(0.5, 0.4, 0.3, 0.2, 0.1)
#' K <- c(0.01, 0.02, 0.03, 0.04, 0.05)
#' lst_neighbours <- matrix(c(2, 3, 1, 3, 4, 2, 4, 5, 3, 5, 4, 2, 5, 3, 1), ncol = 3)
#' IDW_interpolation_nodes(2.5, 5.5, nodes_x, nodes_y, vel_x, vel_y, K, 3, lst_neighbours)
IDW_interpolation_nodes <- function(x,y,nodes_x,nodes_y,vel_x,vel_y,K,id.node,lst_neighbours) {
  
  neighbours_node <- lst_neighbours[id.node,]
  
  nodes_x <- nodes_x[neighbours_node]
  nodes_y <- nodes_y[neighbours_node]
  
  vel_x <- vel_x[neighbours_node]
  vel_y <- vel_y[neighbours_node]
  K <- K[neighbours_node]
  
  distances <- sqrt((nodes_x - x)^2 + (nodes_y-y)^2)
  
  weights <- 1/distances^2
  weights <- weights/sum(weights)
  
  int_vel_x <- sum(weights*vel_x)
  int_vel_y <- sum(weights*vel_y)
  int_K <- sum(weights*K)
  
  node_closest <- neighbours_node[which.min(distances)]

  
  minimum_distance <- as.numeric(distances[which.min(distances)])
  
  return(data.frame(vel_x=int_vel_x,
                    vel_y=int_vel_y,
                    int_K=int_K,
                    node.id=node_closest,
                    minimum_distance=minimum_distance))
  
}