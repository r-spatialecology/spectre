#' @title plot_commonness
#' 
#' @description Plot commonness between observed and optimized data
#' 
#' @param x Results object of run_optimization_min_conf()
#' @param target Pairwise matrix of species in common.
#' 
#' @details 
#' Plot a heatmap of commonness between observed data and optimized data. This visual style allows for easier spatial understanding of commonness differences to be ascertained.
#' 
#' @return ggplot

#' @export
plot_commonness <- function(x, target) {
  species_grid <- x$optimized_grid
  commonness_species_grid <- calculate_solution_commonness_rcpp(species_grid)
  
  n_row <- nrow(commonness_species_grid)
  n_col <- ncol(commonness_species_grid)
  
  commonness_difference_grid <- t(commonness_species_grid - target)
  
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