#' Compute Adaptive Time Step at Each Iteration
#'
#' This function calculates an adaptive time step to use in a simulation based on the distance to the next position and the current flow velocity. The adaptive time step helps improve the accuracy and stability of the simulation by adjusting the time step dynamically based on the current flow conditions.
#'
#' @param distance The distance between the current position and the next position in the simulation.
#' @param vel The flow velocity at the current position.
#' @return The computed adaptive time step for the simulation.
#' @export
#'
#' @examples
#' # Calculate the adaptive time step for a simulation
#' distance <- 2.5 # meters
#' vel <- 0.1 # meters per second
#' adaptive_timestep(distance, vel)
adaptive_timestep <- function(distance,vel,vel_min=0.1) {
  
  if (vel == 0) {
    
    # to avoid dt = Inf if vel = 0
    dt <- 1/4*distance/vel_min
    
  } else {
    
    dt <- 1/4*distance/vel
    
  }
  
  
  
  return(dt)
  
}
