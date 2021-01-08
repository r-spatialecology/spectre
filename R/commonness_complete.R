#' @title commonness_complete
#' 
#' @description Calculates the overall proportion of matching cells between a 
#'  commonness matrix calculated from a species presence absence grid and an 
#'  objective commonness matrix.
#' 
#' @param species_grid A species presence absence matrix. The number of columns
#'  should match the total number of sites in the landscape of interest and the 
#'  number of rows should match the maximum number of species in the area. 
#'  Species presences are denoted with a 1 while absences are denoted with a 0. 
#'  Typically, this matrix would be the output of the 
#'  \code{spectre::run_optimization_min_conf()} function.
#' @param target Pairwise matrix of denoting the number of species in common 
#' between every site in the landscape of interest derived using estimates of
#'  alpha (species richness) and beta (Bray-Curtis dissimilarity) biodiversity
#'  estimates of the landscape. This matrix typically would be the objective
#'  function matrix used as an input for the 
#'  \code{spectre::run_optimization_min_conf()} function. 
#' 
#' @details \code{commonness_complete} initially calculates a pairwise 
#'  commonness matrix from the \code{species_grid} by summing the number
#'  of present species in common for each site-site pair. This "solution"
#'  commonness matrix is then compared with the \code{target} matrix and the
#'  total number of cells with value differences recorded and returned as a
#'  proportion of the total number of cells in each pairwise matrix.
#' 
#' 
#' @return An \code{integer} value of the overall proportion of matching cells
#'  between the calculated "solution" commonness matrix and the \code{target}
#'  commonness matrix.

#' @export
commonness_complete <- function(species_grid, target) {
  
  commonness_species_grid <- calculate_solution_commonness_rcpp(species_grid[[1]])
  
  n_row <- nrow(commonness_species_grid)
  n_col <- ncol(commonness_species_grid)
  
  commonness_species_grid[upper.tri(commonness_species_grid, diag = TRUE)] <- NA
  
  commonness_difference_grid <- commonness_species_grid - target
  
  result <- sum(commonness_difference_grid == 0, na.rm = TRUE) / 
    sum(commonness_species_grid >= 0, na.rm = TRUE)
  
  return(result)
}