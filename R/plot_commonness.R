#' @title plot_commonness
#' 
#' @description Plot commonness
#' 
#' @param species_grid Optimized grid using run_optimization.
#' @param target Pairwise matrix of species in common.
#' 
#' @details 
#' Plot a heatmap of commonnes between observed data and optimized data.
#' 
#' @return ggplot
#' @references xxx

#' @export
plot_commonness <- function(species_grid, target) {
  
  commonness_species_grid <- calculate_solution_commonness_rcpp(species_grid)
  
  n_row <- nrow(commonness_species_grid)
  n_col <- ncol(commonness_species_grid)
  
  commonness_species_grid[upper.tri(commonness_species_grid, diag = TRUE)] <- NA
  
  commonness_difference_grid <- commonness_species_grid - target
  
  commonness_difference <- expand.grid(x = seq_len(n_row), 
                                       y = seq_len(n_col))
  
  commonness_difference$value <- commonness_difference_grid[seq_len(n_row * n_col)]
  
  plot <- ggplot2::ggplot(data = commonness_difference) + 
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value)) + 
    ggplot2::scale_fill_gradient2(name = "Differences in \ncommonness", 
                                  low = "#440154FF", 
                                  mid = "white", 
                                  high = "#FDE725FF") + 
    ggplot2::coord_equal() + 
    ggplot2::theme_void()
  
  return(plot)
}