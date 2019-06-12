#' Function calculate the common species pair-wise matrix (i.e. the target) 
#' I have just made this very simply and it is likely this is an excellent place to try and improve efficiency. 
#' FEELING LIKE I AM MISSING OBVIOUS IMPROVEMENTS HERE

#' 
#' @param alpha predicted alpha diverisity.
#' @param beta pairwise matrix of predicted dissimilarity.
#' @return matrix of common species between each site


calculate_commonness <- function(alpha, beta){
  #Create an empty matrix to fill
  commonness_mat <- matrix(NA, ncol = ncol(beta), nrow = nrow(beta))
  
  for(i in 1:ncol(beta)){
    
   for(j in 1:nrow(beta)){
     commonness_mat[j,i] <- ceiling(((1 - beta[j,i])*(alpha[j] + alpha[i])) / 2)
   } 
  }
  
  return(commonness_mat)
}




