#' @title commonness_complete
#' 
#' @description Counts the number of correctly solved sites
#' 
#' @param species_grid Optimized grid using run_optimization.
#' @param target Pairwise matrix of species in common.
#' 
#' @details 
#' 
#' 
#' @return Proportion of correctly solved sites
#' @references xxx

#' @export
commonness_complete <- function(species_grid, target) {
  
  commonness_species_grid <- calculate_solution_commonness_rcpp(species_grid[[1]])
  
  n_row <- nrow(commonness_species_grid)
  n_col <- ncol(commonness_species_grid)
  
  commonness_species_grid[upper.tri(commonness_species_grid, diag = TRUE)] <- NA
  
  commonness_difference_grid <- commonness_species_grid - target
  
  result <- sum(commonness_difference_grid == 0, na.rm = TRUE) / 
    sum(commonness_species_grid >= 0, na.rm = TRUE)
  
  return(result)
}