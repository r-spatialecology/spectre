# We thankfully downloaded the Barro Colorado Island (BCI) data (Condit et al., 2019) from:
# https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46. 

# We arbitrarily chose tree stem data set number 8 and only need 
# "treeID", "sp" = species name & "quadrat" = grid cell of tree, number xxyy indicates rows and columns

setwd("~/spectre/BCI_tidy/data")

load("bci.tree8.rdata")
data <- bci.tree8[, names(bci.tree8) %in% c("treeID", "sp", "quadrat")]

Sites <- data$quadrat
Species <- data$sp

d <- data.frame(Sites, Species)

(n_sites <- length(unique(Sites))) # 1251 sites 
(n_species <- length(unique(Species))) # 328 species 

# generate a siteXspecies matrix

res <- matrix(nrow = n_species, ncol = n_sites, data = 0)

for (SITES in 1:n_sites){
  print(paste0("Site number: ", SITES))
  for (SPECIES in 1:n_species){
    
    if (sum(which(Sites == unique(Sites)[SITES]) %in% which(Species == unique(Species)[SPECIES])) > 0 ){
      # species SPECIES is present at site SITE, thus 
      res[SPECIES, SITES] <- 1
    }
  }
}

saveRDS(res, file.path(paste0("BCI_tree8_MSP.rds")))
