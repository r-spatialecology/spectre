#' Function calculate the common species pair-wise matrix (i.e. the target)

#' 
#' @param alpha predicted alpha diverisity.
#' @param beta pairwise matrix of predicted dissimilarity.
#' @return matrix of common species between each site


calculate_commonness <- function(alpha, beta){
  #Create an empty matrix to fill
  commonness_mat <- matrix(NA, ncol = ncol(beta), nrow = nrow(beta))
  
  for(i in 1:ncol(beta)){
    
    for(j in 1:nrow(beta)){
      commonness_mat[j,i] <- ceiling(((1 - bmo[j,i])*(alpha_list[j] + alpha_list[i])) / 2)
    } 
  }
  
  return(commonness_mat)
}