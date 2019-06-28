#' @title plot_energy
#' 
#' @description Plot energy
#' 
#' @param species_grid Optimized grid using run_optimization.
#' 
#' @details 
#' Plot energy over time
#' 
#' @return ggplot
#' @references xxx

#' @export
plot_energy <- function(species_grid) {
  
  energy_df <- species_grid[[2]]
  
  plot <- ggplot2::ggplot(data = energy_df) + 
    ggplot2::geom_line(ggplot2::aes(x = i, y = energy)) + 
    ggplot2::labs(x = "Iterations", y = "Energy") +
    ggplot2::theme_bw()
  
  return(plot)
}