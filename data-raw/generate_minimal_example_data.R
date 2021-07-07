### Here we generate data for the minimal example and for test_generate_commonness_matrix_from_gdm.R 

# library("gdm")

# create random species composition
set.seed(42) 

nspecies <- 20
nsites <- 15
presence_prob <- 0.3 # probability for each species to be present at each site

get_species_list <- function(nspecies, nsites, presence_prob)
{
  # generate a siteXspecies list with a random set of species presences/absences
  # fill list with random species presences
  m <- matrix(nrow = nspecies, ncol = nsites, data = 0)
  
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

# random species composition
species_list <- get_species_list(nspecies = nspecies, nsites = nsites, 
                                 presence_prob = presence_prob) 

# calculation of Bray-Curtis dissimilarity with the gdm package: 
# bioData is required by the gdm package, but does not affect 
# the observed Bray-Curtis dissimilarity we will use later
bioData <- data.frame(site_id = 1:nsites, x_coords = rep(13, nsites), 
                      y_coords = rep(10, nsites)) 

bioData <- cbind(bioData, t(species_list)) 

predData <- data.frame(site_id = 1:nsites, preds = runif(nsites))

sitepairs <- gdm::formatsitepair(bioData = bioData, bioFormat = 1, abundance = FALSE, 
                                 siteColumn = "site_id",
                                 XColumn = "x_coords", YColumn = "y_coords", 
                                 predData = predData)

alpha_div <- colSums(species_list) 
gdm_result <- gdm::gdm(sitepairs, geo = TRUE) 

minimal_example_data <- list(species_list = species_list, alpha_list = alpha_div, beta_list = gdm_result$observed)

usethis::use_data(minimal_example_data)
