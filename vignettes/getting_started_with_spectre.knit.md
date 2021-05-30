---
title: "Getting started with `spectre`"
author: "C.E. Simpkins, S. Hanss, M. Hesselbarth, M. Spangenberg and J. Salecker"
date: "2021-05-21"
output: rmarkdown::html_vignette
bibliography: bibliography.bib
vignette: >
  %\VignetteIndexEntry{Getting started with `spectre`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`spectre` is an `R` package which easily implements an advanced optimization algorithm capable of predicting regional community composition at fine spatial resolutions using only sparse biological and environmental data. 
The algorithm underlying `spectre` utilizes estimates of $\alpha$-diversity (i.e. species richness) and $\beta$-diversity (i.e. species dissimilarity) to come up with community composition estimates for all patches within a target region. 
The method used in `spectre` is an adapted version of that presented by Mokany et al. [-@Mokany2011].  

---

# Installation
Install the release version from CRAN:


```r
# install.packages("spectre") #Uncomment when package is on CRAN
```

To install the developmental version of `spectre`, use:

```r
# install.packages("devtools")
# devtools::install_github("r-spatialecology/spectre") #Uncomment when repo is public
```

---

# Use case example
This example acts as a minimal working case and uses simple "simulated" data matching the structure of that needed by the relevant functions. This simple example is used to minimize the time and data storage requirements needed to run this vignette.

### Generating input data
The first step in using the `spectre` package is to gather a estimates for $\alpha$-biodiversity and $\beta$-biodiversity (in the form of Bray-Curtis dissimilarity) for the area of interest at the desired outcome resolution. 
In this example we use a random species composition to i) get an $\alpha$-diversity estimate and ii) to generate a Bray-Curtis dissimilarity estimate in the output format of the `gdm` package [@Fitzpatrick2021]


```r
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

gdm_result <- gdm::gdm(sitepairs, geo = TRUE) 
```

### Running the optimization

We use the input estimates ($\alpha$-diversity and Bray-Curtis dissimilarity) to generate a commonness matrix (i.e. species in common between each site by site pair) using the `generate_commonness_matrix_from_gdm()` function. 
This commonness matrix acts as the objective function (i.e. target) for the optimization algorithm.


```r
# Calculate objective_matrix from (modelled) alpha-diversity and Bray-Curtis dissimilarity
alpha_list <- colSums(species_list) # alpha-diversity of random species community

objective_matrix <- spectre::generate_commonness_matrix_from_gdm(
  gdm_predictions = gdm_result$observed, 
  alpha_list = alpha_list)
```

Once the input estimates and objective function have been obtained the optimization algorithm is straightforward in `spectre`, requiring only one function call. 
Note though that the run time for this function may be high especially for large landscapes with high species diversity and if `max_iterations` is high.


```r
res <- spectre::run_optimization_min_conf(
  alpha_list = alpha_list,
  total_gamma = nspecies,
  target = objective_matrix,
  max_iterations = 1000) # n iterations
```

```
## 
##  > Optimization finished with lowest absolute error = 11 (highest absolute error was: 99 improved by: 88)
```

### Result analysis
`spectre` incorporates functions to allow for easy calculation of certain error metrics, namely the mean absolute commonness error ($MAE_c$) and the relative commonness error ($\% RCE$). $MAE_c$ is the mean of the absolute difference between the solved *solution matrix* and the *objective function*, whereas $\% RCE$ is the $MAE_c$ over the absolute commonness from the *objective function* represented as a percentage.


```r
error_c <- spectre::calc_commonness_error(x = res, objective_matrix = objective_matrix)
```

The *objective function* had a mean commonness of 1.75. 
The mean absolute error between the *objective function* and the solved *solution matrix* was 0.1. 
The solution matrix had an relative commonness error (RCE) of 6%.

These results can be visualized in two ways using functions built into the package. First, one can plot the error of the solved *solution matrix* over time. Second, the commonness error between the final solved *solution matrix* and the *objective function* for each patch can be plotted


```r
# With an increasing number of iterations, the solution matrix improved
spectre::plot_error(x = res)
```

<img src="C:/Users/Craig/Dropbox/Professional/Academic/Gottingen/Manuscripts/spectre/vignettes/getting_started_with_spectre_files/figure-html/unnamed-chunk-7-1.png" width="50%" />

```r
# Plot commonness error between objective function and solution matrix
spectre::plot_commonness(x = res, target = objective_matrix)
```

<img src="C:/Users/Craig/Dropbox/Professional/Academic/Gottingen/Manuscripts/spectre/vignettes/getting_started_with_spectre_files/figure-html/unnamed-chunk-7-2.png" width="50%" />
---

# References