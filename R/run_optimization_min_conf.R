#' @title run_optimization_min_conf
#' 
#' @description Generate an optimized estimate of community composition 
#'  (species presences and absences) for every site in the study area.
#' 
#' @param alpha_list \code{Matrix} of predicted alpha diversity (species richness) in 
#'  each cell.
#' @param total_gamma Total number of species present throughout the entire
#'  landscape.
#' @param target Pairwise matrix of species in common between each site by site 
#'  pair. Only the upper triangle of the matrix is actually needed.
#' @param fixed_species Fixed partial solution with species that are considered 
#'  as given. Those species are not going to be changed during optimization.
#' @param partial_solution An initial \code{matrix} of species presences and 
#'  absences for each site in the landscape. The total number of presences must
#'  match the estimated species richness of each site.
#' @param max_iterations The maximum number of iterations that the optimization
#'  algorithm may run through before stopping.
#' @param seed Seed for random number generator. Seed must be a positive integer value.
#'   \code{seed = NA} means that a random integer is used as seed. 
#' @param verbose If \code{TRUE} (default), a progress report is printed during
#'  the optimization run. 
#' @param interruptible Allow a run to be interrupted before completion.
#' 
#' @details \code{run_optimization_min_conf} is the core function of the 
#'  \code{spectre} package. The underlying algorithm of this function is
#'  adapted from Mokany et al. (2011). A pairwise commonness matrix (having the 
#'  same structure as the \code{target} matrix) is calculated from the 
#'  \code{partial_solution} matrix and the value difference with the 
#'  \code{target} determined. If a difference is present and depending on the 
#'  set stopping criteria the algorithm continues. A random site in the 
#'  presence/absence matrix is selected, and a random presence record at this 
#'  site replaced with an absence. Every absence in the selected site is then 
#'  individually flipped to a presence and the value difference with the 
#'  objective recorded. The presence record which resulted in the lowest value 
#'  difference (minimum conflict) is retained. This cycle continues, with a 
#'  random site selected every iteration, until the pairwise commonness and 
#'  objective matrices match or the algorithm runs beyond the 
#'  \code{max_iterations}.
#' 
#' @return A species presence/absence \code{matrix} of the study landscape.
#' 
#' @references Mokany, K., Harwood, T.D., Overton, J.M., Barker, G.M., &
#'  Ferrier, S. (2011). Combining \eqn{\alpha} and \eqn{\beta} diversity models to fill
#'  gaps in our knowledge of biodiversity. Ecology Letters, 14(10), 1043-1051.

#' @export
run_optimization_min_conf <- function(alpha_list, 
                                      total_gamma, 
                                      target,
                                      max_iterations,
                                      partial_solution = NULL,
                                      fixed_species = NULL,
                                      seed = NA,
                                      verbose = TRUE,
                                      interruptible = TRUE) {
  
  if(is.null(fixed_species)) {
    fixed_species <- matrix()
  }
  
  if (is.null(partial_solution)) {
    partial_solution <- matrix()
  }
  
  if (is.na(seed)) {
    seed <- sample(0:2147483647, 1) # equal to .Machine$integer.max which creates some issue on Windows machines
  }
  
  # we need the upper triangle of the target matrix, only
  # the lower triangle including the diagonal is set NA in the optimizer.
  result = optimizer_min_conf(alpha_list = alpha_list, 
                              total_gamma = total_gamma, 
                              target = target, 
                              max_iterations = max_iterations,
                              partial_solution = partial_solution,
                              fixed_species = fixed_species,
                              seed = seed, 
                              verbose = verbose,
                              interruptible = interruptible)
  
  return(result)
}
