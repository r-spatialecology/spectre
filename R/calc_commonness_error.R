#' @title calc_commonness_error
#' 
#' @description Calculate commonness error
#' 
#' @param x Results object from run_optimization_min_conf.
#' @param objective_matrix Matrix from (modeled) alpha-diversity and Bray-Curtis dissimilarity
#' 
#' @details 
#' Calculate mean absolute commonness error (MAE_c) and relative commonness error in percentage (RCE).
#' 
#' @return vector
#' 
#' @export
calc_commonness_error <- function(x, objective_matrix) {
  
  solution_matrix <- calculate_solution_commonness_rcpp(x$optimized_grid)
  
  # mean absolute commonness error
  MAE_c <- mean(abs(solution_matrix - objective_matrix), na.rm = TRUE) 
  
  # relative commonness error [%]
  RCE <- MAE_c / mean(abs(objective_matrix), na.rm = TRUE) * 100 
  
  result <- c(MAE_c = MAE_c, RCE = RCE)
  
  return(result)
}