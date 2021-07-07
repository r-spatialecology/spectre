# Tests  
# 1. that generate_commonness_matrix_from_gdm() correctly calculates a commonness matrix given siteXsite values for 
# Bray-Curtis dissimilarity and per site richness
# 2. that siteXsite order used in the gdm package is the same as assumed by generate_commonness_matrix_from_gdm(). 
# Since siteXsite order in upcoming versions of the gdm package could change, we test whether orders match. 

# We created a random species composition (15 sites, gamma diversity = 20) and calculated $\alpha$-diversity and 
# Bray-Curtis dissimilarity to be used as input for the `spectre` algorithm 
# (please see `data-raw/generate_minimal_example_data.R` for details). 

# calculate commonness matrix directly from the random species list
obj_commonness <- spectre:::calculate_solution_commonness_rcpp(minimal_example_data$species_list)

# Re-calculate commonness from observed siteXsite Bray-Curtis dissimilarity values and per site richness 
rec_commonness <- spectre:::generate_commonness_matrix_from_gdm(gdm_predictions = minimal_example_data$beta_list, 
                                                                alpha_list = minimal_example_data$alpha_list)
expect_equal(obj_commonness, rec_commonness)


