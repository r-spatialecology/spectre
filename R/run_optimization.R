#' @title run_optimization
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target Pairwise matrix of species in common.
#' @param max_runs Max number of loops before stopping.
#' @param energy_threshold Optimization stops if energy threshold is reached.
#' @param patience Number of runs with no improvement before stopping.
#' @param verbose It TRUE, progress report is printed
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

run_optimization <- function(alpha_list, total_gamma, target, 
                             max_runs, energy_threshold, patience, 
                             verbose = TRUE) {
  
  # get dimensions of matrix
  n_row <- total_gamma
  
  n_col <- length(alpha_list)
  
  # generate initial solution
  current_solution <- matrix(data = 0, 
                             nrow = n_row, 
                             ncol = n_col)
  
  # loop changes the alpha number of species in each site to present (i.e. 1)
  for (n in seq_len(n_col)) {
    
    change_locations <- rcpp_sample(x = seq_len(n_row),
                                    size = alpha_list[n])
    
    current_solution[change_locations, n] <- 1
  }

  # calculate the site x site commonness for the current solution
  solution_commonness <- calculate_solution_commonness_rcpp(current_solution)
  
  solution_commonness[upper.tri(solution_commonness, diag = TRUE)] <- NA
  
  # calculate the difference between target and current solution
  energy <- abs(sum(solution_commonness - target, na.rm = TRUE))
  
  # init patience counter
  unchanged_steps <- 0
  
  # create random col/row ids
  random_col <- rcpp_sample(x = seq_len(n_col),
                            size = max_runs, replace = TRUE)
  
  random_row_1 <- rcpp_sample(x = seq_len(n_row), 
                            size = max_runs, replace = TRUE)
  
  random_row_2 <- rcpp_sample(x = seq_len(n_row), 
                              size = max_runs, replace = TRUE)
  
  # for loop not longer than max_runs
  for (i in seq_len(max_runs)) {
    
    # create a new modified site x species grid
    new_solution <- current_solution

    # get random ids
    current_col <- random_col[i]
    
    current_row_1 <- random_row_1[i]
    
    current_row_2 <- random_row_2[i]

    # change values
    new_solution[current_row_1, current_col] <- current_solution[current_row_2, 
                                                                 current_col]
    
    new_solution[current_row_2, current_col] <- current_solution[current_row_1, 
                                                                 current_col]

    # calculate the site x site commonness for the new solution
    solution_commonness <- calculate_solution_commonness_site_rcpp(new_solution, 
                                                                   solution_commonness, 
                                                                   current_col)

    # calculate the difference between target and new solution
    energy_new <- abs(sum(solution_commonness - target, na.rm = TRUE))

    # check if energy decreased
    if (energy_new < energy) {
      
      # keep new solution
      current_solution <- new_solution
      
      # keep new energy
      energy <- energy_new
      
      # reset patience counter
      unchanged_steps <- 0
    } 
    
    # increase patience counter
    else {
      unchanged_steps <- unchanged_steps + 1
    }
    
    # exit loop if enery threshold or patience counter max is reached
    if (energy <= energy_threshold || unchanged_steps > patience) {
      break
    }
    
    # print progress
    if (verbose) {
      message("\r> Progress: max_runs: ", i, "/", max_runs,
              " || energy = ", energy, "\t\t",
              appendLF = FALSE)
    }
  }
  
  # write result in new line if progress was printed
  if (verbose) {
    message("\r")
    message(paste("> Optimization finished with an energy =", energy))
  }

  return(current_solution)
}