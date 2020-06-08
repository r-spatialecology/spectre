testthat::context("Test min_conf2")
library(dplyr)

testdata <- tibble("1" = c(0,0,0,1,0,0),
                   "2" = c(0,0,1,0,0,0),
                   "3" = c(1,0,0,0,1,0),
                   "4" = c(0,0,1,0,0,1),
                   "5" = c(1,1,1,1,1,0))
alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

res_sim1 <- spectre:::optimizer_min_conf2(alpha_list_test, 
                                          total_gamma_test, 
                                          target_matrix_test,
                                          matrix() ,200, 0.0,
                                          verbose = FALSE)
suppressWarnings(
  resultdata <- res_sim1$optimized_grid %>% as_tibble()
)
alpha_list_result <- resultdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_result <- resultdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result <- spectre:::calculate_solution_commonness_rcpp(res_sim1$optimized_grid)

testthat::expect_equal(alpha_list_test, alpha_list_result)
testthat::expect_equal(total_gamma_test, total_gamma_result)
testthat::expect_equal(target_matrix_test, target_matrix_result)


