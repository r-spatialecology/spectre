#' @title generate_data_virtualspecies
#' 
#' @description Create multiple rasters of virtual species distribution based on suitability maps
#' 
#' @param ncol size of the landscape (columns)
#' @param nrow size of the landscape (rows)
#' @param corr_within strength of spatial autocorrelation within each species
#' @param corr_among defines strength of spatial autocorrelation among species
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
#' Thus, if corr_among = 0, each species will build its own community, resulting in low correlation among species.
#' If corr_among = 1, the same suitability map will be used for each species, resulting in very high correlation among species.
#' 
#' 2. Correlation within species:
#' We also want to investigate the effect of spatial clustering within each species.
#' Thus, we use the gaussian field landscape of NLMR. We can use the correlation parameter to define the strength of spatial autocorrelation.
#' When corr_within = 0, the landscape suitability is distributed more or less random, wheres spatial clustering is present for higher values.
#' 
#' The function will finally report a rasterstack with species distribution (presence/absence) of each species.
#' 
#' @return rasterstack
#' @examples 
#' \dontrun{
#' ncol <- 50
#' nrow <- 50
#' gamma <- 100
#' corr_within <- 0 # low values = little spatial autocorrelation
#' corr_among <- 0.1 # 10% of gamma -> 10 communities
#' beta <- 0.75  # lower beta = wider species distributions
#' 
#' ## Generate species:
#' spp <- generate_data_virtualspecies(ncol, nrow, corr_within, corr_among, gamma, beta)
#' ## Investigate distributions:
#' raster::plot(spp)
#' ## Calculate species richness:
#' alpha <- sum(spp)
#' }
#' @export

generate_data_virtualspecies <- function(ncol,
                                         nrow,
                                         corr_within,
                                         corr_among,
                                         gamma,
                                         beta) {
  
  ## Calculate number of communities:
  communities <- ifelse(corr_among == 0, gamma, round(gamma / (gamma * corr_among)))
  
  ## Generate a suitability map for each community:
  suitability <- purrr::map(seq(communities), function(x){
    NLMR::nlm_gaussianfield(ncol=ncol, nrow=nrow, mag_var = corr_within)}) 
  suitability <- raster::stack(suitability)
  
  ## Generate species distributions:
  spp <- purrr::map(seq(gamma), function(x){
    ## Select community (choose one random suitability distribution)
    id <- sample(raster::nlayers(suitability), 1)
    spp_suitable <- suitability[[id]]
    ## Generate virtual species distribution:
    spp_x <- suppressMessages(virtualspecies::convertToPA(spp_suitable, plot=FALSE, beta=beta))
    spp_x <- spp_x$pa.raster
    return(spp_x)
  })
  spp <- raster::stack(spp)
  
  return(spp)
}