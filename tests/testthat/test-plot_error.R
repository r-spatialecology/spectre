context("plot_error")

library(dplyr)

testdata <- data.frame("1" = c(0,0,0,1,0,0),
                       "2" = c(0,0,1,0,0,0),
                       "3" = c(1,0,0,0,1,0),
                       "4" = c(0,0,1,0,0,1),
                       "5" = c(1,1,1,1,1,0))

alpha_list_test <- testdata %>% summarise_all(sum) %>% as.numeric()
total_gamma_test <- testdata %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- testdata %>% as.matrix() %>% spectre:::calculate_solution_commonness_rcpp()

res_sim <- run_optimization_min_conf(alpha_list = alpha_list_test, 
                                     total_gamma = total_gamma_test, 
                                     target = target_matrix_test, 
                                     max_iterations = 200, 
                                     verbose = FALSE)

test_that("plot_error returns ggplot", {
  
  plot_e <- plot_error(x = res_sim)

  expect_is(object = plot_e, class = "ggplot")
})

test_that("plot_commonness returns ggplot", {
  
  plot_c <- plot_commonness(x = res_sim, target = target_matrix_test)
  
  expect_is(object = plot_c, class = "ggplot")
})

