#' @title generate_data_raster_solution
#' 
#' @description Create a solution matrix from a species richness raster landscape
#' 
#' @param landscape any raster that represents species richness as cell values (e.g. derived from generate_data_raster_rescale())
#' @param gamma Total (estimated) species in the system.
#' 
#' @details 
#' This function takes a raster landscape which contains species richness as cell values. 
#' Then, an empty matrix is created with one row for each species (gamma) and one column for each site (ncol * nrow of landscape).
#' Then, on each cell, the raster value defines the number of present species IDs (rows). 
#' The selection of rows to be present, is chosen completely at random.
#' 
#' Such a solution matrix can be used to calculate an optimization target matrix using the calculate_solution_commonness_rcpp() function (see extended examples).
#' 
#' @return Matrix (species X sites)
#' @examples 
#' \dontrun{
#' # 1. Generate random landscape with NLMR:
#' landscape <- NLMR::nlm_random(ncol=50, nrow=50) 
#' # 2. Rescale to derive species richness values:
#' gamma <- 500
#' alpha <- generate_data_raster_rescale(landscape, gamma, 0.1, 0.4)
#' landscapetools::show_landscape(alpha)
#' # 3. Calculate solution matrix:
#' solution <- generate_data_raster_solution(alpha, gamma)
#' 
#' 
#' ### EXTENDED EXAMPLES:
#' #
#' # Here are some examples for a complete data generation workflow, using NLMR landscape generation functions:
#' 
#' # Define landscape properties:
#' landscape_size <- 50
#' gamma <- 100
#' alpha.min <- 0.1 # minimum species richness (fraction of gamma)
#' alpha.max <- 0.5 # maximum species richness (fraction of gamma)
#' 
#' ## Example 1: Random uniform landscape
#' landscape <- NLMR::nlm_random(ncol=landscape_size, nrow=landscape_size) 
#' alpha <- generate_data_raster_rescale(landscape, gamma, alpha.min, alpha.max)
#' landscapetools::show_landscape(alpha)
#' solution <- generate_data_raster_solution(alpha, gamma)
#' target <- calculate_solution_commonness_rcpp(solution)
#' 
#' ## Example 2: Gaussian field landscape (spatial autocorrelation)
#' landscape <- NLMR::nlm_gaussianfield(ncol=landscape_size, nrow=landscape_size)
#' alpha <- generate_data_raster_rescale(landscape, gamma, alpha.min, alpha.max)
#' landscapetools::show_landscape(alpha)
#' solution <- generate_data_raster_solution(alpha, gamma)
#' target <- calculate_solution_commonness_rcpp(solution)
#' 
#' ## Example 3: Merge different landscape types:
#' landscape_high_autocorrelation <- NLMR::nlm_edgegradient(ncol = landscape_size, nrow = landscape_size, direction = 80)
#' landscape_low_autocorrelation <- NLMR::nlm_fbm(ncol = landscape_size, nrow = landscape_size, fract_dim = 0.5)
#' landscape <- landscapetools::util_merge(landscape_low_autocorrelation, landscape_high_autocorrelation)
#' alpha <- generate_data_raster_rescale(landscape, gamma, alpha.min, alpha.max)
#' landscapetools::show_landscape(alpha)
#' solution <- generate_data_raster_solution(alpha, gamma)
#' target <- calculate_solution_commonness_rcpp(solution)
#' }
#' @export

generate_data_raster_solution <- function(landscape, gamma)
{
  ## Convert to site by species matrix:
  current_solution <- matrix(data = 0, 
                             nrow = gamma, 
                             ncol = raster::ncell(landscape))
  
  ## Select random species in each site:
  for(i in 1:ncol(current_solution))
  {
    species.ids <- unique(as.integer(sample(1:gamma, raster::values(landscape)[i])))
    current_solution[species.ids, i] = 1
  }
  
  return(current_solution)
}
