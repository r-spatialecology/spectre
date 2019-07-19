#' Function calculate the common species pair-wise matrix (i.e. the target)

#' 
#' @param alpha predicted alpha diverisity.
#' @param beta pairwise matrix of predicted dissimilarity.
#' @return matrix of common species between each site


calculate_target <- function(alpha, beta){
  #Create an empty matrix to fill
  target_mat <- matrix(NA, ncol = ncol(beta), nrow = nrow(beta))
  
  for(i in 1:ncol(beta)){
    
    for(j in 1:nrow(beta)){
      target_mat[j,i] <- ceiling(((1 - beta[j,i])*(alpha_list[j] + alpha_list[i])) / 2)
    } 
  }
  
  target_mat[upper.tri(target_mat, diag = TRUE)] <- NA
  
  return(target_mat)
}