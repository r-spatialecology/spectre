#' @title run_optimization
#' 
#' @description xxx
#' 
#' @param alpha_list Matrix of predicted alpha diversity in each cell.
#' @param total_gamma Total (estimated) species in the system.
#' @param target Pairwise matrix of species in common.
#' @param max_runs Max number of loops before stopping.
#' @param annealing Probability to keep modified matrix even if energy increased.
#' @param energy_threshold Optimization stops if energy threshold is reached. This is set as a value between 0 and 1 determining the proportion of error accepted
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
run_optimization <- function(alpha_list, 
                             total_gamma, 
                             target, 
                             max_runs,
                             annealing = 0.1,
                             energy_threshold, 
                             patience, 
                             verbose = TRUE) {
  
  # get dimensions of matrix
  n_row <- total_gamma
  n_col <- length(alpha_list)
  
  # generate initial solution
  current_solution <- matrix(data = 0, 
                             nrow = n_row, 
                             ncol = n_col)
  
  # generate dataframe for i and energy
  energy_df <- data.frame(i = rep(NA, max_runs), 
                          energy = rep(NA, max_runs))
  
  # loop changes the alpha number of species in each site to present (i.e. 1)
  for (n in seq_len(n_col)) {
    
    change_locations <- rcpp_sample(x = seq_len(n_row),
                                    size = alpha_list[n])
    
    current_solution[change_locations, n] <- 1
  }

  # calculate the site x site commonness for the current solution
  solution_commonness <- calculate_solution_commonness_rcpp(current_solution)
  
  # Not necessary?
  # solution_commonness[upper.tri(solution_commonness, diag = TRUE)] <- NA
  
  # calculate the difference between target and current solution
  energy <- sum(solution_commonness != target, na.rm = TRUE) / sum(!is.na(target), na.rm = TRUE)
  
  # init patience counter
  unchanged_steps <- 0
  
  # # create random col/row ids
  random_col <- rcpp_sample(x = seq_len(n_col),
                            size = max_runs, replace = TRUE)

  # random_row_1 <- rcpp_sample(x = seq_len(n_row), 
  #                           size = max_runs, replace = TRUE)
  # 
  # random_row_2 <- rcpp_sample(x = seq_len(n_row), 
  #                             size = max_runs, replace = TRUE)
  
  # create random number for annealing probability
  if (annealing != 0) {
    
    annealing_random <- stats::runif(n = max_runs, min = 0, max = 1)
  } 
  
  else {
    
    annealing_random <- rep(0, max_runs)
  }
  
  # for loop not longer than max_runs
  for (i in seq_len(max_runs)) {
    
    # create a new modified site x species grid
    new_solution <- current_solution

    # get random col id
    current_col <- random_col[i]
    
    # get random row 1
    current_row_1 <- rcpp_sample(x = seq_len(n_row),
                                 size = 1)
    
    # value of row 1
    value_1 <- current_solution[current_row_1, 
                                current_col]
    
    # sample row where value != value_1
    current_row_2 <- rcpp_sample(x = which(current_solution[, current_col] != value_1), 
                                 size = 1)
    
    # value_2 opposite to value 1
    value_2 <- ifelse(test = value_1 == 1, yes = 0, no = 1)
    
    # change values
    new_solution[current_row_1, current_col] <- value_2
    
    new_solution[current_row_2, current_col] <- value_1

    # calculate the site x site commonness for the new solution
    solution_commonness <- calculate_solution_commonness_site_rcpp(new_solution, 
                                                                   solution_commonness, 
                                                                   current_col)

    # calculate the difference between target and new solution
    energy_new <- sum(solution_commonness != target, na.rm = TRUE) / sum(!is.na(target), na.rm = TRUE)

    # check if energy decreased
    if (energy_new < energy || annealing_random[i] < annealing) {
      
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
    
    # save energy in df
    energy_df[i, ] <- c(i, energy)
    
    # print progress
    if (verbose) {
      message("\r> Progress: max_runs: ", i, "/", max_runs,
              " || energy = ", round(energy, 5), "\t\t\t",
              appendLF = FALSE)
    }
  }
  
  # write result in new line if progress was printed
  if (verbose) {
    message("\r")
    message(paste("> Optimization finished with an energy =", round(energy, 5)))
  }
  
  result <- list(current_solution, 
                 energy_df)
  
  names(result) <- c("optimized_grid", "energy")

  return(result)
}