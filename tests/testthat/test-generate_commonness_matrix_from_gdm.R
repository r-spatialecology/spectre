# tests whether the generate_commonness_matrix_from_gdm function works
library("gdm")
# Create random siteXspecies data
nspecies <- 15
nsites <- 10
presence_prob <- 0.3 # probability for each species to be present at each site

get_objective_matrix <- function(nspecies, nsites, presence_prob)
{
  # generate a random matrix with n species and n sites
  # presence_prob is the probability for each species to appear on a site
  # It's only a quick and dirty approach here :-)
  m <- matrix(nrow = nsites, ncol = nspecies, data = 0)
  
  # fill matrix with random species presences
  
  for (row in 1:ncol(m)){
    for (col in 1:nrow(m)){
      if (runif(1) < presence_prob){
        m[col, row] <- 1
      }
    }
  }
  mode(m) <- "integer"
  return(m)
}


(obj_matrix <- get_objective_matrix(nspecies = nspecies, nsites = nsites, presence_prob = presence_prob))
obj_matrix[1:2, ] <- 1

alpha_list <- rowSums(obj_matrix)

m <- as.data.frame(obj_matrix)
names(obj_matrix) <- paste0("spec_", 1:nspecies)

(obj_commonness <- spectre:::calculate_solution_commonness_rcpp(t(obj_matrix)))

bioData <- data.frame(site_id = 1:nsites, x_coords = rep(13, nsites), y_coords = rep(10, nsites))

bioData <- cbind(bioData, obj_matrix)

predData <- data.frame(site_id = 1:nsites, preds = runif(nsites))


sitepairs <- gdm::formatsitepair(bioData = bioData, bioFormat = 1, abundance = FALSE, 
                    siteColumn = "site_id",
                    XColumn="x_coords", YColumn="y_coords", 
                    predData = predData)
haha <- gdm::gdm(sitepairs, geo = TRUE)
haha$observed

spectre:::generate_commonness_matrix_from_gdm(haha$observed, alpha_list = alpha_list)
