#' Simulate Trajectory with Nearest Neighbor Interpolation
#'
#' This function simulates the trajectory of a larva in a flow field based on given
#' initial positions and flow data. The larva's trajectory is computed by updating
#' its position iteratively according to the flow velocities at each step, until
#' certain conditions are met (e.g., reaching a "safe zone" or a boundary).
#'
#' Nearest neighbor interpolation is used to find the flow velocities and kinematic
#' properties (e.g., turbulence) at each iteration. The interpolation is performed using
#' the provided `lst_neighbours`, which represents the matrix of nearest neighbors
#' for each node in the flow field.
#'
#' @param flow_df A data frame containing flow data, including `x`, `y`, `vel_x`, `vel_y`,
#'                `K`, `depth`, `boundary`, and other relevant columns.
#' @param pos0 A list or data frame representing the initial position of the larva,
#'             containing `x`, `y`, and other relevant columns.
#' @param n.it An integer specifying the number of iterations for simulating the trajectory.
#' @param up A numeric value representing the vertical velocity of the larva.
#' @param vel_coef A numeric coefficient used to correct the velocity based on position.
#' @param lst_neighbours A matrix representing the nearest neighbors for each node in the flow
#'                       field. This matrix is used for nearest neighbor interpolation.
#' @return A list containing the simulated trajectory of the larva (`traj_larvae`),
#'         information about the drift destination (`traj_info`), and the final position
#'         of the larva (`posf`).
#' @import dplyr
#' @export
#'
#' @examples
#' # Define flow data and initial position
#' flow_df <- data.frame(
#'   x = c(1, 2, 3, 4, 5),
#'   y = c(2, 4, 6, 8, 10),
#'   vel_x = c(0.1, 0.2, 0.3, 0.4, 0.5),
#'   vel_y = c(0.5, 0.4, 0.3, 0.2, 0.1),
#'   K = c(0.01, 0.02, 0.03, 0.04, 0.05),
#'   depth = c(1, 2, 3, 4, 5),
#'   boundary = c("dry", "wet", "wet", "outflow", "dry")
#' )
#' pos0 <- list(x = 1, y = 2, id = 1)
#' n.it <- 100
#' up <- 0.1
#' lst_neighbours <- matrix(c(2, 3, 1, 3, 4, 2, 4, 5, 3, 5, 4, 2, 5, 3, 1), ncol = 3)
#' simulate_trajectory2(flow_df, pos0, n.it, up, lst_neighbours)
simulate_trajectory <- function(flow_df,pos0,n.it,up,vel_coef=1,lst_neighbours) {
  
  library(dplyr)
  
  # correct velocity based on coefficient
  flow_df <- flow_df %>%
    dplyr::mutate(vel_x = vel_x*vel_coef,
                  vel_y = vel_y*vel_coef,
                  K = K*vel_coef) # assuming that also effect of turbulence is proportional to position in the vertical (and scales as velocity)
  
  # initialize vectors
  nodes_k <- numeric(length=n.it)
  
  # initialize positions with NA, so it's easier to identify early interruptions of drift (when for loop interrupted with break)
  
  dt_k <- rep(NA,n.it)
  
  ux_k <- rep(NA,n.it)
  uy_k <- rep(NA,n.it)
  
  vel_k <- rep(NA,n.it)
  depth_k <- rep(NA,n.it)
  
  upx_k <- numeric(length=n.it)
  upy_k <- numeric(length=n.it)
  
  x_k <- numeric(length=n.it)
  y_k <- numeric(length=n.it)
  
  K_k <- numeric(length=n.it)
  
  # first step - take data from closest node
  x_k[1] <- pos0$x
  y_k[1] <- pos0$y
  
  # find id of nearest node
  distance <- sqrt((flow_df$x - x_k[1])^2 + (flow_df$y - y_k[1])^2)
  node.id <- which.min(distance)[1]
  dist_min <- as.numeric(distance[node.id])
  
  # interpolate
  interp <- IDW_interpolation_nodes(x_k[1],y_k[1],flow_df$x,flow_df$y,flow_df$vel_x,flow_df$vel_y,flow_df$K,node.id,lst_neighbours)
  
  nodes_k[1] <- node.id
  
  ux_k[1] <- flow_df$vel_x[node.id]
  uy_k[1] <- flow_df$vel_y[node.id]
  
  ux_k[1] <- interp$vel_x
  uy_k[1] <- interp$vel_y
  
  depth_k[1] <- flow_df$depth[node.id]
  vel_k[1] <- sqrt(ux_k[1]^2 + uy_k[1]^2)
  
  dt_k[1] <- adaptive_timestep(dist_min,vel_k[1])
  
  K_k[1] <- interp$int_K
  
  drift_destination <- 'in flow' # in case larvae doesn't reach any of the boundaries
  outside <- FALSE
  outside_n <- 0
  
  # start iterations
  for (i in 2:n.it) {
    
    # compute new larvae  position
    pos <- cpp_update_position_inflow(x_k[i-1],y_k[i-1],ux_k[i-1],uy_k[i-1],up,K_k[i-1],dt_k[i-1],drift_destination)
    
    x_k[i] <- pos$x
    y_k[i] <- pos$y
    
    # interpolate flow velocities at larvae position
    interp <- IDW_interpolation_nodes(x_k[i],y_k[i],flow_df$x,flow_df$y,flow_df$vel_x,flow_df$vel_y,flow_df$K,node.id,lst_neighbours)
    
    dist_min <- interp$minimum_distance
    node.id <- interp$node.id
    
    nodes_k[i] <- node.id
    
    ux_k[i] <- interp$vel_x
    uy_k[i] <- interp$vel_y
    
    K_k[i] <- interp$int_K
    
    depth_k[i] <- flow_df$depth[node.id]
    vel_k[i] <- sqrt(ux_k[i]^2 + uy_k[i]^2)
    
    dt_k[i] <- adaptive_timestep(dist_min,vel_k[i])
    
    # check if larvae in "nursery" (an area with velocity < 10 cm/s), outside boundary, or outflow
    
    n.it_max <- i
    
    drift_destination <- 'in flow'
    
    if (flow_df$boundary[node.id] == 'dry') {
      drift_destination <- 'outside boundary'
      
      # save info if particle ever been outside boundary
      outside <- TRUE
      outside_n <- outside_n + 1
      
      if (outside_n > 1000) {
        break
      }
      
      # reset position
      x_k[i] <- x_k[i-1]
      y_k[i] <- y_k[i-1]
      
    } else if (vel_k[i] <= 0.1) {
      drift_destination <- 'nursery'
      break
    } else if (flow_df$boundary[node.id] == 'outflow') {
      drift_destination <- 'outflow'
      break
    }
    
  }
  
  traj_larvae <- data.frame(data.frame(x=x_k,
                                       y=y_k,
                                       ux=ux_k,
                                       uy=uy_k,
                                       vel=vel_k,
                                       depth=depth_k,
                                       dt=dt_k,
                                       node=nodes_k))
  traj_larvae <- traj_larvae[1:n.it_max,] # select only up to maximum computed iteration
  
  return(list(traj_larvae=traj_larvae,
              traj_info=data.frame(drift_destination=drift_destination,
                                   outside=outside,
                                   n.it_max=n.it_max,
                                   id=pos0$id),
              posf=data.frame(x=traj_larvae$x[n.it_max],
                              y=traj_larvae$y[n.it_max])))
}