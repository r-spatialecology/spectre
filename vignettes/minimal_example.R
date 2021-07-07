## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("spectre")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("r-spatialecology/spectre")

## -----------------------------------------------------------------------------
library("spectre")

## -----------------------------------------------------------------------------
# load "observed" alpha-, beta- and gamma-diversity values of the random species composition
alpha_list <- minimal_example_data$alpha_list # richness
beta_list <- minimal_example_data$beta_list # Bray-Curtis dissimilarity
total_gamma <- dim(minimal_example_data$species_list)[1] # 20 species

## ---- message = FALSE, warning = FALSE----------------------------------------
library("spectre")

objective_matrix <- spectre::generate_commonness_matrix_from_gdm(gdm_predictions = beta_list, 
                                                                 alpha_list = alpha_list)

# Solve composition 
res <- spectre::run_optimization_min_conf(alpha_list = alpha_list,
                                          total_gamma = total_gamma,
                                          target = objective_matrix,
                                          max_iterations = 1000, # n iterations
                                          seed = 123) # use a random seed for reproducibility

## ---- message = FALSE---------------------------------------------------------
error_c <- spectre::calc_commonness_error(x = res, objective_matrix = objective_matrix)

## ---- include = TRUE, out.width="50%", fig.align="center"---------------------
# With an increasing number of iterations, the solution matrix improved
spectre::plot_error(x = res)

# Plot commonness error between objective matrix and solution matrix
spectre::plot_commonness(x = res, target = objective_matrix)

