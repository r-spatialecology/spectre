#' @title generate_data_virtualspecies
#' 
#' @description Create multiple rasters of virtual species distribution based on suitability maps
#' 
#' @param landscape_size size (sidelength) of the landscape
#' @param correlation strength of spatial autocorrelation within each species
#' @param community_size size of species communities, defined as proportion of gamma (defines strength of spatial autocorrelation between species: large communities = strong correlation)
#' @param gamma Total (estimated) species in the landscape
#' @param beta Species suitability treshhold
#' 
#' @details 
#' 
#' The function generates virtual species distributions, using the virtualspecies package.
#' These species are distributed according to suitability maps that are created with the NLMR package first.
#' In order to provide both, control over spatial autocorrelation within and between species the following approach is applied:
#' 
#' 1. Correlation among species:
#' We assume, that a community has a certain habitat suitability. Each community however, may consist of multiple species.
#' Thus, when creating sitability maps, we do not neccessarily create one for each species, but one for each community.
#' Thus, if community_size = 0, each species will build its own community, resulting in low correlation among species.
#' If community_size = 1, the same suitability map will be used for each species, resulting in very high correlation among species.
#' 
#' 2. Correlation within species:
#' We also want to investigate the effect of spatial clustering within each species.
#' Thus, we use the gaussian field landscape of NLMR. We can use the correlation parameter to define the strength of spatial autocorrelation.
#' When correlation = 0, the landscape suitability is distributed more or less random, wheres spatial clustering is present for higher values.
#' 
#' The function will finally report a rasterstack with species distribution (presence/absence) of each species.
#' 
#' @return rasterstack
#' @examples 
#' 
#' landscape_size <- 50
#' gamma <- 100
#' correlation <- 0 # low values = little spatial autocorrelation
#' community_size <- 0.1 # 10% of gamma -> 10 communities
#' beta <- 0.75  # lower beta = wider species distributions
#' 
#' ## Generate species:
#' spp <- generate_data_virtualspecies(landscape_size, correlation, community_size, gamma, beta)
#' ## Investigate distributions:
#' raster::plot(spp)
#' ## Calculate species richness:
#' alpha <- sum(spp)
#' 
#' @export

generate_data_virtualspecies <- function(landscape_size,
                                         correlation,
                                         community_size,
                                         gamma,
                                         beta) {
  
  ## Calculate number of communities:
  communities <- ifelse(community_size == 0, gamma, round(gamma / (gamma * community_size)))
  
  ## Generate a suitability map for each community:
  suitability <- purrr::map(seq(communities), function(x){
    NLMR::nlm_gaussianfield(ncol=landscape_size, nrow=landscape_size, mag_var = correlation)}) 
  suitability <- raster::stack(suitability)
  
  ## Generate species distributions:
  spp <- purrr::map(seq(gamma), function(x){
    ## Select community (choose one random suitability distribution)
    id <- sample(raster::nlayers(suitability), 1)
    spp_suitable <- suitability[[id]]
    ## Generate virtual species distribution:
    spp_x <- virtualspecies::convertToPA(spp_suitable, plot=FALSE, beta=beta)
    spp_x <- spp_x$pa.raster
    return(spp_x)
  })
  spp <- raster::stack(spp)
  
  return(spp)
}