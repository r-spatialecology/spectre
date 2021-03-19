# Try simulated annealing with spectre
# slim version to evalate performace of T_0 and temp_stop = 

# Analogue to verbose, see how parameters are handled in the package... 


out_dir <- "./results/SimAnnB01"
getwd()
dir.create(file.path(out_dir))

library("spectre")
library("foreach")
library("tidyr")
doParallel::registerDoParallel(12)

# Parameters
# ===========
replicate <- 1:10 # 1:15
nspecies <- 50
nsites <- 100
max_iterations <- c(c(200) * 1000) # > 150,000 not needed 
acceptance_prob <- seq(0, 0.005, .0005) # probability to accept ~ 100 of all intermediate solutions, given temperature = T_0

parameters <- tidyr::crossing(replicate,
                              nspecies,
                              nsites,
                              max_iterations,
                              acceptance_prob)

parameters$seed <- 1:dim(parameters)[1] 

sim_fun <- function(siminputrow, parameters, writeRDS, verbose)
{
  # Read and set parameters
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) {
    p <- parameters[10, ] # for quick testing
  }
  replicate <- p$replicate
  
  
  seed <- p$seed
  set.seed(seed)
  
  
  
  
  
  ### 
  time_before <- Sys.time()
  
  nspecies <- 30
  nsites <- 20
  max_iterations <- 50000
  T_0 <- 50*1000 
  temp_stop <- 30000 
  
  # Generate objective and a starting solution
  obj_list <- spectreSimAnn::get_objective_matrix(nspecies = nspecies,
                                                  nsites = nsites,
                                                  presence_prob = 0.2)
  
  target <- spectre:::calculate_solution_commonness_rcpp(obj_list)
  alpha_list <- colSums(obj_list)
  
  
  time_before <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_list, 
                                                     total_gamma = nspecies, 
                                                     target = target, 
                                                     max_iterations = max_iterations,
                                                     seed = seed,
                                                     verbose = TRUE,
                                                     interruptible = TRUE,
                                                     T_0 = 1000) 
  iteration <- res_min_conf$error$i
  error <- res_min_conf$error$error
  
  plot(error ~ iteration, pch = 20, main = paste0("T_0= ", T_0, "; temp_stop= ", temp_stop), ylim = c(0, max(error)), xlim = c(0, 2000))
  
  time_after <- Sys.time()
  (calc_time <- time_after - time_before)
  
  solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
  MAE_c <- mean(abs(solution_commonness - target), na.rm = TRUE) 
  (RCE <- MAE_c / mean(abs(target), na.rm = TRUE) * 100 )
  
  r_1 <- tibble::tibble(replicate = replicate,
                        n_sites = nsites,
                        n_species = nspecies,
                        MAE_c = MAE_c,
                        RCE = RCE,
                        interations = max_iterations,
                        time = calc_time,
                        acceptance_prob = acceptance_prob)
  
  saveRDS(r_1, file = paste0(out_dir, "/res", siminputrow, "f.rds")) 
}

print(paste0(dim(parameters)[1], " replicates"))

foreach(REPLICATE = 1:dim(parameters)[1], .export = c("p"), .packages = c("spectre")) %dopar% {
  print(paste0("Replicate: ", REPLICATE))
  
  sim_fun(REPLICATE,
          parameters,
          TRUE,
          TRUE)
  return(1)
} 

