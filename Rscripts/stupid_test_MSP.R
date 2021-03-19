
nspecies <- 50
nsites <- 30
max_iterations <- 20000
T_0_factor <- 0.0001
# temp_stop <- 30000 

# Generate objective and a starting solution
obj_list <- spectreSimAnn::get_objective_matrix(nspecies = nspecies,
                                                nsites = nsites,
                                                presence_prob = 0.2)

target <- spectre:::calculate_solution_commonness_rcpp(obj_list)
alpha_list <- colSums(obj_list)

solution_list  <- spectreSimAnn::get_solution_matrix(nspecies = nspecies,
                                                     nsites = nsites,
                                                     richness_vector = alpha_list)

(T_0 <- spectreSimAnn::calculate_T_0(objective =  target,
                                     solution_list = solution_list,
                                     n_samples = 2000,
                                     acceptance_prob = 0.8))
T_0 <- T_0 * T_0_factor

res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_list, 
                                                   total_gamma = nspecies, 
                                                   target = target, 
                                                   max_iterations = max_iterations,
                                                   seed = 1,
                                                   verbose = TRUE,
                                                   interruptible = TRUE,
                                                   T_0 = T_0)  
iteration <- res_min_conf$error$i
error <- res_min_conf$error$error

plot(error ~ iteration, pch = 20, main = paste0("T_0_factor= ", round(T_0_factor, 4)), ylim = c(0, max(error)) ) # , xlim = c(0, 2000))
