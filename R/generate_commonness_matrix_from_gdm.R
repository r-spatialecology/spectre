#' @title generate_commonness_matrix_from_gdm
#' 
#' @description Generates the best estimate objective (siteXsite commonness matrix) from predicted richness and predicted Bray-Curtis dissimilarity. 
#' 
#' 
#' @details The function expects predicted per site richness as a vector, and Bray-Curtis dissimilarity as provided by the gdm-package (Fitzpatrick et al. 2020) output (as a list).
#' @return 
#' @examples example
#' @export
### Generate commonness matrix from sorensen gdm-package output 

generate_commonness_matrix_from_gdm <- function(gdm_predictions, alpha_list){
  # takes predictions from the gdm-package (in list form) and 
  # generates the best estimate objective matrix
  
  n_sites <- length(alpha_list)
  
  if( !(length(gdm_predictions) == (n_sites^2 - n_sites))){ # quick check
    print(paste0("Probably number of sites in estimated richness and estimated gdm predictions do not match! Please check inputs"))
  }
  
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
  
  target_commonness <- matrix(nrow = n_sites, ncol = n_sites, data = NA)
  for (site in 1:n_sites){
    for (other_site in 1:n_sites){
      if (other_site <= site){
        target_commonness[other_site, site] <- NA
      } else {
        beta <- target_sorensen[other_site, site]
        target_commonness[other_site, site] <- round((1 - beta) * (alpha_list[site] + alpha_list[other_site]) / 2 )
      }
      
    }
  }
  
  return(target_commonness) 
}
