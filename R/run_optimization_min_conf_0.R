#' @title run_optimization_min_conf_0
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target Pairwise matrix of species in common.
#' @param fixed_species Fixed partial solution with species that are considered as given. Those species are not going to be changed during optimization.
#' @param partial_solution Partial or complete initial solution as a start for the optimization.
#' @param max_iterations Max number of loops before stopping.
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
run_optimization_min_conf_0 <- function(alpha_list, 
                                        total_gamma, 
                                        target,
                                        fixed_species = NULL,
                                        partial_solution = NULL,
                                        max_iterations,
                                        energy_threshold,
                                        seed = 0,
                                        verbose = TRUE) {
  
  if(is.null(fixed_species)) {
    fixed_species <- matrix()
  }
  if (is.null(partial_solution)) {
    partial_solution <- matrix()
  }
  result = optimizer_min_conf0(alpha_list = alpha_list, 
                               total_gamma = total_gamma, 
                               target = target, 
                               fixed_species,
                               partial_solution = partial_solution,
                               max_iterations = max_iterations,
                               energy_threshold = energy_threshold,
                               seed = seed, 
                               verbose = verbose)
  
  return(result)
}
