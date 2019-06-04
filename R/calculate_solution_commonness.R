#' Function calculates the commonness matrix for site-species solutions 
#' Uses equation: Cij = sum(total_species(i) * total_species(j)) 
#' AGAIN I THINK THIS IS A PLACE WHERE WE CAN IMPROVE EFFICIENCY
 
#' @param solution_matrix predicted species-site matrix.
#' @return matrix of common species between each site

calculate_solution_commonness <- function(solution_matrix){
  ## Create an empty matrix to fill
  commonness_s_mat <- matrix(NA, ncol = ncol(solution_matrix), nrow = ncol(solution_matrix))
  
  for(i in 1:ncol(commonness_s_mat)){
    
    for(j in 1:nrow(commonness_s_mat)){
      commonness_s_mat[j,i] <- sum((solution_matrix[,j] + solution_matrix[,i]) > 1)
    } 
  }
  
  ## Remove the lower triangle of the matrix to match target
  commonness_s_mat[upper.tri(commonness_s_mat, diag = TRUE)] <- NA
  
  return(commonness_s_mat)
}