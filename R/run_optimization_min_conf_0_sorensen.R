#' @title run_optimization_min_conf_0_sorensen
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target_sorensen Pairwise matrix of sorensen dissimilarity.
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
#' @return Best matrix of Sorensen dissimilarity between each site.
#' @references xxx

#' @export
run_optimization_min_conf_0_sorensen <- function(alpha_list, 
                                                 total_gamma, 
                                                 target_sorensen, # MSP
                                                 fixed_species = NULL,
                                                 partial_solution = NULL,
                                                 max_iterations,
                                                 patience = 2000,
                                                 energy_threshold,
                                                 seed = 0,
                                                 verbose = TRUE,
                                                 norm = "sum",
                                                 p = 1) {
  
  if(is.null(fixed_species)) {
    fixed_species <- matrix()
  }
  if (is.null(partial_solution)) {
    partial_solution <- matrix()
  }
  
  result = optimizer_min_conf0_sorensen(alpha_list = alpha_list, 
                                        total_gamma = total_gamma, 
                                        target_sorensen = target_sorensen,# MSP
                                        fixed_species = fixed_species,
                                        partial_solution = partial_solution,
                                        max_iterations = max_iterations,
                                        patience = patience,
                                        energy_threshold = energy_threshold,
                                        seed = seed, 
                                        verbose = verbose,
                                        norm = norm,
                                        p = p)
  
  return(result)
}
