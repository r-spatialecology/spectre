context("calc_commonness_error")

library(dplyr)

testdata <- data.frame("1" = c(0,0,0,1,0,0),
                       "2" = c(0,0,1,0,0,0),
                       "3" = c(1,0,0,0,1,0),
                       "4" = c(0,0,1,0,0,1),
                       "5" = c(1,1,1,1,1,0))

alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

target_matrix_test[3, 4] <- 2
target_matrix_test[3, 5] <- 0

res_sim <- run_optimization_min_conf(alpha_list = alpha_list_test, 
                                     total_gamma = total_gamma_test, 
                                     target = target_matrix_test, 
                                     max_iterations = 200, 
                                     verbose = FALSE)

test_that("calc_commonness_error works", {
  
  error_c <- calc_commonness_error(x = res_sim, objective_matrix = target_matrix_test)
  
  expect_true(all(error_c != 0))
  expect_is(object = error_c, class = "numeric")
  expect_length(object = error_c, n = 2)
  expect_named(object = error_c, expected = c("MAE_c", "RCE"))
})
