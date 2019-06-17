#' @title run_optimization
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target pairwise matrix of species in common.
#' @param cycles number of loops before stopping.
#' @param required_D maximum D required to stop.
#' @param patience number of loops with no improvement before stopping.
#' 
#' @details 
#' This is the function which runs the optimization algorithm.
#' There are a number of steps involved
#' 1. generate initial solution site x species matrix (only needs to be run once per landscape) 
#' 3. alter solution site x species matrix and select best site x site matrix based on Dsol
#' 4. Run optimization algorithm until some stop requirement is met 
#' 
#' @return Best matrix of common species between each site.
#' @references xxx
#' @export
#' 

run_optimization <- function(alpha_list, total_gamma, target, cycles, required_D, patience){
  
  ## 1) Generate initial solution
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  current_solution <- matrix(0, nrow = total_gamma, ncol = length(alpha_list))

  #Loop changes the alpha number of species in each site to present (i.e. 1)
  for(n in 1:ncol(current_solution)){
    #print(alpha_list[n])
    change_locations <- sample(nrow(current_solution),alpha_list[n])
    current_solution[change_locations, n] <- 1
  }

  #Calculate the site x site commonness for the current solution
  solution_commonness <- calculate_solution_commonness_rcpp(current_solution)
  solution_commonness[upper.tri(solution_commonness, diag = TRUE)] <- NA
  
  #Calculate the difference between target and current solution
  D <- abs(sum(solution_commonness - as.numeric(target), na.rm = TRUE))

  #While loop to change matrix and check new D till best solution found
  #END REQUIREMENTS: 1) D < x; 2) n > y; 3) delta D not changing over set period
  n_loops <- 0
  unchanged_steps <- 0
  while (n_loops < cycles && D > required_D && unchanged_steps < patience) {
    #print(D)
    #Create a new modified site x species grid
    new_solution <- current_solution

    random_col <- sample(seq_len(ncol(new_solution)), size = 1)
    random_row1 <- sample(seq_len(nrow(new_solution)), size = 1)
    random_row2 <- sample(seq_len(nrow(new_solution)), size = 1)

    new_solution[random_row1, random_col] <- current_solution[random_row2, random_col]
    new_solution[random_row2, random_col] <- current_solution[random_row1, random_col]

    #Calculate the site x site commonness for the new solution
    solution_commonness <- calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col)

    #Calculate the difference between target and new solution
    D_new <- abs(sum(solution_commonness - as.numeric(target), na.rm = TRUE))

    if(D_new < D){
      current_solution <- new_solution
      D <- D_new
      unchanged_steps <- 0
    } else {
      unchanged_steps <- unchanged_steps + 1
    }
    n_loops <-  n_loops + 1
  }

  print(paste("finished with D of", D))

  return(current_solution)
}