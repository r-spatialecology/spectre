library(dplyr)

testdata <- matrix(c(0,0,0,1,0,
                     0,0,1,0,0,
                     1,0,0,0,1,
                     0,0,1,0,0,
                     1,1,1,1,1), nrow = 5, ncol = 5)
species_prop <- c(2/5, 1/5, 3/5, 2/5, 2/5)

alpha_list_test <- testdata %>% as_tibble() %>% summarise_all(funs(sum)) %>% slice(1) %>% unlist(., use.names = FALSE)
total_gamma_test <- testdata %>% as_tibble() %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
target_matrix_test <- spectre:::calculate_solution_commonness_rcpp(testdata)

species_grid <- run_optimization(alpha_list = alpha_list_test,
                                 total_gamma = total_gamma_test,
                                 target = target_matrix_test,
                                 max_runs = 2000,
                                 energy_threshold = 0.0,
                                 patience = 5000)

test <- spectre:::mh_optimizer(alpha_list_test, total_gamma_test, target_matrix_test, species_prop, 20000, 0.0)
#test2 <- spectre:::mh_optimizer_neutral(alpha_list_test, total_gamma_test, target_matrix_test, species_prop, 20000)
mat <- test[[1]]
target_matrix_test2 <- spectre:::calculate_solution_commonness_rcpp(mat)
spectre:::plot_commonness(test, target_matrix_test)
spectre:::plot_energy(test)
sum(abs((target_matrix_test - target_matrix_test2)), na.rm = TRUE) /
  sum((target_matrix_test2), na.rm = TRUE)
spectre:::calc_energy(target_matrix_test, target_matrix_test2)



# Data gerneration for unit tests
species_grid <- run_optimization(alpha_list = spectre::alpha_list,
                                 total_gamma = spectre::estimated_gamma,
                                 target = spectre::target_matrix,
                                 max_runs = 20000,
                                 energy_threshold = 0.1,
                                 annealing = 0.1,
                                 patience = 5000)
species_grid_test <- species_grid$optimized_grid

# Unit test, if alpha diversity does not change
alpha_list_test <- species_grid_test %>% as_tibble() %>% summarise_all(funs(sum)) %>% slice(1) %>% unlist(., use.names=FALSE)
identical(spectre::alpha_list, alpha_list_test)

# Unit test, if gamma diversity does not change
total_gamma_test <- species_grid_test %>% as_tibble() %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()
identical(139L, total_gamma_test)

# Unit test, if solution is good enough
target_matrix_test <- spectre:::calculate_solution_commonness_rcpp(species_grid_test)

species_grid <- run_optimization(alpha_list = alpha_list_test,
                                 total_gamma = total_gamma_test,
                                 target = target_matrix_test,
                                 max_runs = 150000,
                                 energy_threshold = .1,
                                 patience = 500)

species_grid_test2 <- species_grid$optimized_grid
target_matrix_test2 <- spectre:::calculate_solution_commonness_rcpp(species_grid_test2)
energy <- spectre:::calc_energy(target_matrix_test2, target_matrix_test)

energy < 0.25  # or whatever

test <- spectre:::mh_optimizer(spectre::alpha_list, 139L, target_matrix, 20000, base_probability_jump = 0.1)
test2 <- spectre:::mh_optimizer_neutral(spectre::alpha_list, 139L, spectre::target_matrix, 20000)

plot_commonness(species_grid = test,
                target = target_matrix)
plot_energy(species_grid = test2)
plot_energy(species_grid)
plot_energy(species_grid = test2)
plot_energy(species_grid = test)
