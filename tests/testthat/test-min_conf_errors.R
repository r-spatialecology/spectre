# tests whether an error is thrown when the partial solution and fixed species parameters do not fit
context("min_conf_errors")

library(dplyr)

testdata <- data.frame("1" = c(0,0,0,1,0,0),
                       "2" = c(0,0,1,0,0,0),
                       "3" = c(1,0,0,0,1,0),
                       "4" = c(0,0,1,0,0,1),
                       "5" = c(1,1,1,1,1,0))

partial_solution_fixed <- data.frame("1" = c(0,0,0,0,0,0),
                                     "2" = c(0,0,0,0,0,0),
                                     "3" = c(0,0,0,0,0,0),
                                     "4" = c(0,0,0,0,0,0),
                                     "5" = c(1,1,1,1,1,1)) %>% as.matrix()

partial_solution <- data.frame("1" = c(0,0,0,0,0,0),
                               "2" = c(0,0,0,0,0,0),
                               "3" = c(0,0,0,0,0,0),
                               "4" = c(0,0,0,0,0,0)) %>% as.matrix()

fixed_species <- data.frame("1" = c(0,0,0,0,0,0),
                            "2" = c(0,0,1,0,0,0),
                            "3" = c(1,0,0,0,1,0),
                            "4" = c(0,0,1,0,0,1),
                            "5" = c(1,1,1,1,1,1),
                            "6" = c(1,0,0,0,1,0)) %>% as.matrix()

alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% calculate_solution_commonness_rcpp()

expect_error(
  run_optimization_min_conf(alpha_list = alpha_list_test, 
                            total_gamma = total_gamma_test, 
                            target = target_matrix_test,
                            partial_solution = partial_solution,
                            max_iterations = 0, 
                            verbose = FALSE),
  regexp = "The size of the partial_solution vector does not match n_sites * gamma_div. partial_solution ignored.", 
  fixed = TRUE)

expect_error(
  run_optimization_min_conf(alpha_list = alpha_list_test, 
                            total_gamma = total_gamma_test, 
                            target = target_matrix_test,
                            partial_solution = partial_solution_fixed,
                            fixed_species = fixed_species,
                            max_iterations = 0, 
                            verbose = FALSE),
  "The size of the fixed_species vector does not match n_sites * gamma_div. fixed_species ignored.", 
  fixed = TRUE)
