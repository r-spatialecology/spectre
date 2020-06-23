#' @title generate_data_raster_rescale
#' 
#' @description Utility function for raster based data generator to rescale raster values as species richness
#' 
#' @param landscape any raster, for example generated with the NLMR package
#' @param gamma Total (estimated) species in the system.
#' @param alpha.min Minimum nunmber of species in a cell (defined as a fraction of gamma, e.g. 0.1 represents 10% of gamma)
#' @param alpha.max Maximum nunmber of species in a cell (defined as a fraction of gamma, e.g. 0.5 represents 50% of gamma)
#' 
#' @details 
#' This function rescales any raster landscape with numerical values to the desired biodiversity indicators.
#' It uses the following formula for rescaling:
#' NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
#' 
#' 
#' @return Rescaled raster landscape
#' @examples 
#' \dontrun{
#' # 1. Generate random landscape with NLMR:
#' landscape <- NLMR::nlm_random(ncol=50, nrow=50) 
#' # 2. Rescale to derive species richness values:
#' gamma <- 500
#' alpha <- generate_data_raster_rescale(landscape, gamma, 0.1, 0.4)
#' landscapetools::show_landscape(alpha)
#' }
#' @export

generate_data_raster_rescale <- function(landscape, 
                                       gamma, 
                                       alpha.min, 
                                       alpha.max)
{
  # NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  landscape_rescaled <- round((((landscape - raster::cellStats(landscape, "min")) * (alpha.max - alpha.min)) / 
                         (raster::cellStats(landscape, "max") - raster::cellStats(landscape, "min")) + alpha.min) * gamma)
  return(landscape_rescaled)
}