#' @title calculate_rce
#' 
#' @description Calculate the relative commonness error
#' 
#' @param solution_list Site x species matrix of the spectre-generated solution.
#' @param target Pairwise matrix of species in common of the objective.
#' 
#' @details 
#' 
#' 
#' @return rce 
#' @references xxx

#' @export
calculate_rce <- function(solution_list, target)
{
  solution <- spectre:::calculate_solution_commonness_rcpp(solution_list)
  error <- mean(abs(solution - target), na.rm = TRUE)
  mean_target_commonness <- mean(target, na.rm = TRUE)
  RCE <- error / mean_target_commonness * 100  
  return(RCE)
}