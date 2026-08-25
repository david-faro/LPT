#' Compute Adaptive Time Step
#'
#' Computes an adaptive time step based on the distance to the next position
#' and the current flow velocity. If the flow velocity is zero, a minimum
#' velocity is used to avoid an infinite time step.
#'
#' @param distance Numeric. Distance to the next position, in meters.
#' @param vel Numeric. Flow velocity at the current position, in meters per
#'   second.
#' @param vel_min Numeric. Minimum flow velocity, in meters per second, used
#'   to calculate the time step when `vel` is zero. Defaults to `0.1`.
#'
#' @return A numeric value representing the adaptive time step, in seconds.
#'
#' @export
#'
#' @examples
#' # Calculate the adaptive time step
#' adaptive_timestep(distance = 2.5, vel = 0.1)
#'
#' # Use a custom minimum velocity when velocity is zero
#' adaptive_timestep(distance = 2.5, vel = 0, vel_min = 0.05)
adaptive_timestep <- function(distance, vel, vel_min = 0.1) {

  if (vel == 0) {

    # Avoid dt = Inf when vel = 0
    dt <- 1 / 4 * distance / vel_min

  } else {

    dt <- 1 / 4 * distance / vel

  }

  return(dt)
}
