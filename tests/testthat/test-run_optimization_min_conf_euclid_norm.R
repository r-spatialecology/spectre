testthat::context("Test min_conf euclid norm")

library(dplyr)

testdata <- tibble("1" = c(0,0,0,1,0,0),
                   "2" = c(0,0,1,0,0,0),
                   "3" = c(1,0,0,0,1,0),
                   "4" = c(0,0,1,0,0,1),
                   "5" = c(1,1,1,1,1,0))
alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

res_sim0 <- run_optimization_min_conf_0(alpha_list = alpha_list_test, 
                                        total_gamma = total_gamma_test, 
                                        target = target_matrix_test, 
                                        max_iterations = 200, 
                                        energy_threshold = 0.0,
                                        verbose = FALSE,
                                        norm = "euclid")

suppressWarnings(
  resultdata0 <- res_sim0$optimized_grid %>% as_tibble()
)

alpha_list_result0 <- resultdata0 %>% summarise_all(sum) %>% as.numeric()
total_gamma_result0 <- resultdata0 %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result0 <- calculate_solution_commonness_rcpp(res_sim0$optimized_grid)

testthat::expect_equal(alpha_list_test, alpha_list_result0)
testthat::expect_equal(total_gamma_test, total_gamma_result0)
testthat::expect_equal(target_matrix_test, target_matrix_result0)