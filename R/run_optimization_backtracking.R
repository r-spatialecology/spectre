#' @title run_optimization
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target Pairwise matrix of species in common.
#' @param max_runs Max number of loops before stopping.
#' @param energy_threshold Optimization stops if energy threshold is reached. This is set as a value between 0 and 1 determining the proportion of error accepted
#' @param seed Seed for random number generator. seed = 0 means that a time stamp is used as seed. 
#' @param verbose It TRUE, progress report is printed
#' 
#' @details 
#' This is the function which runs the optimization algorithm.
#' 
#' @return Best matrix of common species between each site.
#' @references xxx

#' @export
run_optimization_backtracking <- function(alpha_list, 
                                      total_gamma, 
                                      target, 
                                      max_runs,
                                      energy_threshold,
                                      seed = 0,
                                      verbose = TRUE) {
  
  result = optimizer_backtracking(alpha_list = alpha_list, 
                                  total_gamma = total_gamma, 
                                  target = target, 
                                  max_iterations = max_runs,
                                  seed = seed, 
                                  verbose = verbose)
  
  return(result)
}