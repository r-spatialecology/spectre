#' Find best solution
#'
#' @param alpha_list 
#' @param total_gamma 
#' @param target 
#' @param max_iterations 
#' @param partial_solution 
#' @param fixed_species 
#' @param seed 
#' @param max_repetitions 
#' @param stop_at_optimum 
#'
#' @return
#' @export
#'
#' @examples
find_solutions_min_conf <- function(alpha_list,
                                    total_gamma,
                                    target,
                                    max_iterations,
                                    partial_solution = NULL,
                                    fixed_species = NULL,
                                    seed = NA,
                                    max_repetitions,
                                    stop_at_optimum = TRUE) {
  
  best_error <- .Machine$integer.max
  error <- integer()
  results <- list()
  
  if (is.na(seed)) {
    seed <- rep(NA, times = max_repetitions)
  } else {
    seed <- seed:(seed+max_repetitions-1)
  }
  
  for (i in 1:length(seed)) {
    res <- run_optimization_min_conf(alpha_list = alpha_list,
                                     total_gamma = total_gamma,
                                     target = target,
                                     max_iterations = max_iterations,
                                     partial_solution = partial_solution,
                                     fixed_species = fixed_species,
                                     seed = seed[i], verbose = FALSE, interruptible = FALSE)
    
    results[[i]] <- res
    error <- c(error, min(res$error$error))
    if (min(res$error$error) < best_error[i]) {
      best_error <- c(best_error, min(res$error$error))
      if (stop_at_optimum) (
        if (min(res$error$error) == 0) {
          break
        }
      )
    } else {
      best_error <- c(best_error,  best_error[i])
    }
  }
  
  return(list(seed = seed[1:length(error[-1])], best_error = error[-1], error = error, results = results))
}
