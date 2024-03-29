# tests whether the partial solution and fixed species parameters work
context("fixed_species")

library(dplyr)

testdata <- data.frame("1" = c(0,0,0,1,0,0),
                       "2" = c(0,0,1,0,0,0),
                       "3" = c(1,0,0,0,1,0),
                       "4" = c(0,0,1,0,0,1),
                       "5" = c(1,1,1,1,1,0))

partial_solution <- data.frame("1" = c(0,0,0,0,0,0),
                               "2" = c(0,0,0,0,0,0),
                               "3" = c(0,0,0,0,0,0),
                               "4" = c(0,0,0,0,0,0),
                               "5" = c(1,1,1,1,1,1)) %>% as.matrix()

fixed_species <- data.frame("1" = c(1,1,1,1,1,1),
                            "2" = c(0,0,1,0,0,0),
                            "3" = c(1,0,0,0,1,0),
                            "4" = c(0,0,1,0,0,1),
                            "5" = c(1,1,1,1,1,1)) %>% as.matrix()

alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

res_sim <- run_optimization_min_conf(alpha_list = alpha_list_test, 
                                     total_gamma = total_gamma_test, 
                                     target = target_matrix_test,
                                     fixed_species = fixed_species,
                                     partial_solution = partial_solution,
                                     max_iterations = 200, 
                                     verbose = FALSE)

resultdata <- res_sim$optimized_grid %>% as.data.frame()

alpha_list_result <- resultdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_result <- resultdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_result <- spectre:::calculate_solution_commonness_rcpp(res_sim$optimized_grid)

resultdata <- resultdata %>% as.matrix()

expect_equal(alpha_list_result[5], 6L)
expect_equal(alpha_list_result[1], 0L)
expect_equal(total_gamma_test, total_gamma_result)
expect_true(resultdata[1,3] == 0L)
expect_true(resultdata[3,2] == 0L)
expect_true(resultdata[3,4] == 0L)
expect_true(resultdata[5,3] == 0L)
expect_true(resultdata[6,4] == 0L)
expect_setequal(resultdata[1:6,5], rep(1L, times = 6))
