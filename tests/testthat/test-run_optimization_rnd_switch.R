testthat::context("Test random switch")
library(dplyr)

testdata <- matrix(c(0,0,0,1,0,0,
                     0,0,1,0,0,0,
                     1,0,0,0,1,0,
                     0,0,1,0,0,1,
                     1,1,1,1,1,0), nrow = 6, ncol = 5)
alpha_list_test <- testdata %>% as_tibble() %>% summarise_all(funs(sum)) %>% slice(1) %>% unlist(., use.names = FALSE)
total_gamma_test <- testdata %>% as_tibble() %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- spectre:::calculate_solution_commonness_rcpp(testdata)

res_sim1 <- spectre::run_optimization_rnd_switch(alpha_list = alpha_list_test, 
                                                 total_gamma = total_gamma_test, 
                                                 target = target_matrix_test, 
                                                 max_runs = 2000, 
                                                 energy_threshold = 0.0,
                                                 patience = 500,
                                                 verbose = FALSE)
resultdata <- res_sim1$optimized_grid

alpha_list_result <- resultdata %>% as_tibble() %>% summarise_all(funs(sum)) %>% slice(1) %>% unlist(., use.names = FALSE)
total_gamma_result <- resultdata %>% as_tibble() %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result <- spectre:::calculate_solution_commonness_rcpp(resultdata)

testthat::expect_equal(alpha_list_test, alpha_list_result)
testthat::expect_equal(total_gamma_test, total_gamma_result)
testthat::expect_equal(target_matrix_test, target_matrix_result)


