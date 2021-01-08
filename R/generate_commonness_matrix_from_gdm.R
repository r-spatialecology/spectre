#' @title generate_commonness_matrix_from_gdm
#' 
#' @description Creates a pairwise site by site commonness matrix from 
#'  estimates of species richness and Bray-Curtis dissimilarity.
#' 
#' @details \code{generate_commonness_matrix_from_gdm} uses a vector of 
#'  estimated species richness per site and a pairwise matrix of site by site
#'  Bray-Curtis dissimilarity (we recommend using the gdm-package
#'  (Fitzpatrick et al. 2020) to generate this matrix) to produce a matrix of 
#'  the estimated species in common between site pairs (referred to as a
#'  commonness matrix). The commonness between sites is  calculated using
#'  \deqn{C_{ij}=(1-\beta_{ij})(S_{i} + S_{j})/2}
#'  Where \eqn{\beta_{ij}} is the dissimilarity between sites, \eqn{C_{ij}} is
#'  the species in common between sites, and S is the number of species in 
#'  each site. For more details see Mokany et al 2011.
#'
#' @references Mokany, K., Harwood, T.D., Overton, J.M., Barker, G.M., &
#'  Ferrier, S. (2011). Combining \alpha and \beta diversity models to fill
#'  gaps in our knowledge of biodiversity. Ecology Letters, 14(10), 1043-1051. 
#'
#' @return A pairwise site by site \code{matrix} of the number of species in
#'  common between each site pair, with dimensions equal to that of the 
#'  provided dissimilarity matrix.
#'  
#' @param gdm_predictions a square pairwise \code{matrix} of Bray-Curtis 
#'  dissimilarity estimates between site pairs. We recommend using the 
#'  gdm-package (Fitzpatrick et al. 2020) to generate this matrix
#' 
#' @param alpha_list a \code{vector} of species richness for every site in the
#'  study area. The length of this vector must be equivalent to one of the
#'  dimensions of the \code{gdm_predictions}
#' 
#' @export
### Generate commonness matrix from sorensen gdm-package output 

generate_commonness_matrix_from_gdm <- function(gdm_predictions, alpha_list){
  # takes predictions from the gdm-package (in list form) and 
  # generates the best estimate objective matrix
  
  n_sites <- length(alpha_list)
  
  if( !(length(gdm_predictions) == ( (n_sites^2 - n_sites) / 2 ) ) ){ # quick check
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
