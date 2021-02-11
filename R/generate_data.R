### Generate target commonness matrix from randomly initiated site x species matrix  

# generate random  ; loop changes the alpha number of species in each site to present (i.e. 1)

generate_data_simple <- function(total_gamma, n_sites, alpha_list) {
  
  current_solution <- matrix(data = 0, 
                             nrow = total_gamma, 
                             ncol = n_sites)
  
  for (n in seq_len(n_sites)) {
    
    change_locations <- sample(x = seq_len(total_gamma),
                               size = alpha_list[n])
    
    current_solution[change_locations, n] <- 1
  }
  
  target <- calculate_solution_commonness_rcpp(current_solution)
  return(target)
}

get_sample_points <- function(successful_iterations, n_sample_points){
  sample_points <- c(1, seq(round(successful_iterations/n_sample_points), (successful_iterations -1), round(successful_iterations/n_sample_points))                      ,successful_iterations)
    return(sample_points)
}
