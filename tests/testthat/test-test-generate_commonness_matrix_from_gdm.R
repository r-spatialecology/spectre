# Tests  
# 1. that generate_commonness_matrix_from_gdm() correctly calculates a commonness matrix given siteXsite values for 
# Bray-Curtis dissimilarity and per site richness
# 2. that siteXsite order used in the gdm package is the same as assumed by generate_commonness_matrix_from_gdm(). 
# Since siteXsite order in upcoming versions of the gdm package could change, we test whether orders match. 

library("gdm")

set.seed(42) # use seed to ensure the gdm package will converge with a given random species composition and arbitrary set predictors

# Create random siteXspecies data
nspecies <- 50
nsites <- 40
presence_prob <- 0.3 # probability for each species to be present at each site

get_objective_matrix <- function(nspecies, nsites, presence_prob)
{
  # generate a random matrix with n species and n sites
  # fill matrix with random species presences
  
  m <- matrix(nrow = nsites, ncol = nspecies, data = 0)
  
  for (row in 1:ncol(m)) {
    for (col in 1:nrow(m)) {
      if (runif(1) < presence_prob) {
        m[col, row] <- 1
      }
    }
  }
  mode(m) <- "integer"
  return(m)
}

# Create random species composition 
obj_matrix <- get_objective_matrix(nspecies = nspecies, nsites = nsites, presence_prob = presence_prob)

alpha_list <- rowSums(obj_matrix) # richness per site 

# Commonness matrix of the random siteXspecies matrix, used for evaluation later
obj_commonness <- spectre:::calculate_solution_commonness_rcpp( t(obj_matrix) )

# Bray-Curtis dissimilarity is calculated using the gdm package
bioData <- data.frame(site_id = 1:nsites, x_coords = rep(13, nsites), y_coords = rep(10, nsites))
bioData <- cbind(bioData, obj_matrix)

predData <- data.frame(site_id = 1:nsites, preds = runif(nsites))

sitepairs <- gdm::formatsitepair(bioData = bioData, bioFormat = 1, abundance = FALSE, 
                                 siteColumn = "site_id",
                                 XColumn = "x_coords", YColumn = "y_coords", 
                                 predData = predData)
gdm_result <- gdm::gdm(sitepairs, geo = TRUE)

# Re-calculate commonness from observed siteXsite Bray-Curtis values and per site richness 
rec_commonness <- spectre:::generate_commonness_matrix_from_gdm(gdm_result$observed, alpha_list = alpha_list)

expect_equal(obj_commonness, rec_commonness)
