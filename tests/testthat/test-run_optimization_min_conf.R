testthat::context("Test min_conf")

library(dplyr)

testdata <- tibble("1" = c(0,0,0,1,0,0),
                   "2" = c(0,0,1,0,0,0),
                   "3" = c(1,0,0,0,1,0),
                   "4" = c(0,0,1,0,0,1),
                   "5" = c(1,1,1,1,1,0))
alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

res_sim1 <- spectre::run_optimization_min_conf(alpha_list = alpha_list_test, 
                                        total_gamma = total_gamma_test, 
                                        target = target_matrix_test, 
                                        max_runs = 200, 
                                        energy_threshold = 0.0,
                                        verbose = FALSE)
res_sim2 <- spectre::run_optimization_min_conf_1(alpha_list = alpha_list_test, 
                                               total_gamma = total_gamma_test, 
                                               target = target_matrix_test, 
                                               max_runs = 200, 
                                               energy_threshold = 0.0,
                                               verbose = FALSE)
res_sim3 <- spectre::run_optimization_min_conf_2(alpha_list = alpha_list_test, 
                                               total_gamma = total_gamma_test, 
                                               target = target_matrix_test, 
                                               max_runs = 200, 
                                               energy_threshold = 0.0,
                                               verbose = FALSE)
suppressWarnings(
  resultdata <- res_sim1$optimized_grid %>% as_tibble()
)
suppressWarnings(
resultdata1 <- res_sim2$optimized_grid %>% as_tibble()
)
suppressWarnings(
resultdata2 <- res_sim3$optimized_grid %>% as_tibble()
)


alpha_list_result <- resultdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_result <- resultdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result <- spectre:::calculate_solution_commonness_rcpp(res_sim1$optimized_grid)

testthat::expect_equal(alpha_list_test, alpha_list_result)
testthat::expect_equal(total_gamma_test, total_gamma_result)
testthat::expect_equal(target_matrix_test, target_matrix_result)

alpha_list_result1 <- resultdata1 %>% summarise_all(sum) %>% as.numeric()
total_gamma_result1 <- resultdata1 %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result1 <- spectre:::calculate_solution_commonness_rcpp(res_sim2$optimized_grid)

testthat::expect_equal(alpha_list_test, alpha_list_result1)
testthat::expect_equal(total_gamma_test, total_gamma_result1)
testthat::expect_equal(target_matrix_test, target_matrix_result1)

alpha_list_result2 <- resultdata2 %>% summarise_all(sum) %>% as.numeric()
total_gamma_result2 <- resultdata2 %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result2 <- spectre:::calculate_solution_commonness_rcpp(res_sim3$optimized_grid)

testthat::expect_equal(alpha_list_test, alpha_list_result2)
testthat::expect_equal(total_gamma_test, total_gamma_result2)
testthat::expect_equal(target_matrix_test, target_matrix_result2)
