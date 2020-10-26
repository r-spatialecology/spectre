#' @title generate_sorensen_matrix_from_gdm
#' 
#' @description Generate a sorensen siteXsite matrix from gdm-package output
#' 
#' 
#' @details details
#' @return 
#' @examples example
#' @export
### Generate target sorensen matrix from gdm-package-output 

generate_sorensen_matrix_from_gdm <- function(gdm_predictions, n_sites){
  # takes predictions from the gdm-package (in list form) and 
  # generates a target sorensen matrix
  
  target_sorensen <- matrix(nrow = n_sites, ncol = n_sites, data = NA)
  
  INDEX <- 1 
  
  for(i in 1:n_sites){ # 1 von 1:5 
    j <- i + 1
    while(j <= n_sites){
      target_sorensen[j, i] <- gdm_predictions[INDEX]
      INDEX <- INDEX + 1
      j <- j + 1 
    }
  }
  return(target_sorensen) 
}