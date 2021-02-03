# Test sim anneal with dev branch
# the last days I implemented sim anneal in spectre branch sim_ann_3
# in spectreSimAnn, T_0 is calculated as mean(error)-acceptance 

library("foreach")
library("tidyr")
doParallel::registerDoParallel(8)

out_dir <- "./results/002_SA002"
getwd()
dir.create(file.path(out_dir))

### Parameters:
replicate <- 1:10
nspecies <- 100
nsites <- 50
max_iterations <- 100000
presence_prob <- c(0.2, 0.3) # probability for each species to be present at each site
T_0_factor <- seq(0, 0.00025, length.out = 11) # to scale T_0 T_0 < 0.001 was best 

p <- tidyr::crossing(replicate,
                     nspecies,
                     nsites,
                     max_iterations,
                     presence_prob,
                     T_0_factor)

p$seed <- 1:dim(p)[1]

sim_fun <- function(siminputrow, parameters){
  
  # Parameter setting
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) {
    p <- p[1, ] # for quick testing
  }
  replicate <- p$replicate
  nspecies <- p$nspecies
  nsites <- p$nsites
  max_iterations <- p$max_iterations
  presence_prob <- p$presence_prob # probability for each species to be present at each site
  T_0_factor <- p$T_0_factor
  seed <- p$seed
  set.seed(seed)
  
  # Generate objective and a starting solution
  obj_list <- spectreSimAnn::get_objective_matrix(nspecies = nspecies,
                                                  nsites = nsites,
                                                  presence_prob = presence_prob)
  
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
  
  start_time <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_list, 
                                                     total_gamma = nspecies, 
                                                     target = target, 
                                                     max_iterations = max_iterations,
                                                     seed = 1,
                                                     verbose = TRUE,
                                                     interruptible = TRUE,
                                                     T_0 = T_0)  
  stop_time <- Sys.time()
  calc_time <- stop_time - start_time
  
  # With simulated annealing, the minimal error is not neccessarily the final error, thus we save minimal error also
  # in practice, always saving the the "best solution" may not be feasible, bu it could still be interesting... 
  min_error <- min(res_min_conf$error$error)
  (RCE <- spectre:::calculate_rce(res_min_conf$optimized_grid, target))
  
  r_1 <- tibble::tibble(replicate = replicate,
                        nspecies = nspecies,
                        nsites = nsites,
                        RCE = RCE,
                        max_iterations = max_iterations,
                        presence_prob = presence_prob,
                        T_0 = T_0, 
                        min_error = min_error,
                        T_0_factor = T_0_factor,
                        calc_time = calc_time)
  saveRDS(r_1, file = paste0(out_dir, "/comm", siminputrow, ".rds"))
  
}

dim(p)
print(paste0("Start time = ", Sys.time()))

foreach(REPLICATE = 1:dim(p)[1], .export = c("p")) %dopar% {
  
  sim_fun(REPLICATE,
          p)
}
print(paste0("Stopping time = ", Sys.time()))
