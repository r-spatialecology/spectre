#' @title generate_data_virtualspecies_to_solution
#' 
#' @description Transforms rasterstack containing presence/absence maps to a solution matrix
#' 
#' @param species raster stack containing presence/absence maps
#' 
#' @details 
#' 
#' The function transforms a raster stack containing presence/absence maps of species to a solution matrix.
#' This can for example be used for species derive from the virtual species function generate_data_virtualspecies().
#' The solution matrix can then be used to run the spectre optimization algorithms.
#' 
#' @return matrix (species X sites)
#' @examples 
#' \dontrun{
#' ## This example presents the complete workflow of generating virtual species distributions
#' ## over trasnfromation to solution matrix
#' ## to execution of spectre optimization algorithms:' 
#' 
#' 
#' #### Step 1: Generate virtualspecies
#' 
#' # Define parameters:
#' ncol <- 50
#' nrow <- 50
#' gamma <- 100
#' corr_within <- 0 # low values = little spatial autocorrelation
#' corr_among <- 0.1 # 10% of gamma -> 10 communities
#' beta <- 0.75  # lower beta = wider species distributions
#' 
#' ## Generate species:
#' spp <- generate_data_virtualspecies(ncol, nrow, corr_within, corr_among, gamma, beta)
#' alpha <- getValues(sum(spp))
#' 
#' #### Step 2: Transform to solution matrix and calculate target matrix:
#' 
#' solution <- generate_data_virtualspecies_to_solution(spp)
#' target <- calculate_solution_commonness_rcpp(solution)
#' 
#' #### Step 3: Execute optimization
#' 
#' result <- spectre::run_optimization_min_conf(alpha, gamma, target, max_runs = 20000)
#' }
#' @export

generate_data_virtualspecies_to_solution <- function(spp) {
  ## Convert to site by species matrix:
  current_solution <- matrix(data = 0, 
                             nrow = raster::nlayers(spp), 
                             ncol = raster::ncell(spp))
  
  for(i in seq(raster::nlayers(spp)))
  {
    sp_present <- raster::getValues(spp[[i]])
    current_solution[i,] <- sp_present
  }
  return(current_solution)
} 