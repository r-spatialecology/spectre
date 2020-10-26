# Test min_conf_algorithms 

library("spectre")
library("foreach")


### 
sim_fun <- function(siminputrow, parameters, writeRDS, verbose)  # Code mostly stolen from Jan Salecker... 
{
  # Read and set parameters
  p <- parameters[siminputrow, ]
  n_sites <- p$n_sites
  mean_alpha <- p$mean_alpha 
  richness_vector <- p$richness_vector
  richness_sd <- p$richness_sd 
  total_gamma <- p$total_gamma
  replicate <- p$replicate
  max_runs <- p$max_runs
  n_sample_points <- p$n_sample_points
  energy_threshold <- p$energy_threshold 
  tabu <- round(p$tabu_percent * p$n_sites / 100)
  
  # Set random seed
  set.seed(replicate)
  
  # generate target matrix, used for all algorithms
  alpha_list <- round(rnorm(n = n_sites, 
                            mean = mean_alpha, 
                            sd = richness_sd))
  
  # check that richness per site is > 0 and <= gamma diversity. 
  alpha_list[alpha_list < 1] <- 1 # only positive species numbers allowed 
  alpha_list[alpha_list > (total_gamma -2) ] <- (total_gamma -2)  
  
  target_matrix_sim <- spectre:::generate_data_simple(total_gamma = total_gamma, 
                                                     n_sites = n_sites, 
                                                     alpha_list = alpha_list)
  
  
  
  #print("start")
  # Executing original algorithm:
  
  start_time <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf_0(alpha_list = alpha_list, 
                                                     total_gamma = total_gamma, 
                                                     target = target_matrix_sim, 
                                                     max_iterations = max_runs,
                                                     patience = 2500, 
                                                     energy_threshold = energy_threshold,
                                                     # tabu = tabu,
                                                     verbose = verbose) # 
  end_time <- Sys.time()
  time_min_conf <- as.numeric(end_time - start_time)
  
  # only keep a fraction of results (first, sample_points in between , last)
  successful_iterations <- sum(!is.na(res_min_conf$energy[2]))
  # print(paste0("Successful iterations: ", successful_iterations))
  
  if(successful_iterations >= n_sample_points){
    sample_points <- spectre:::get_sample_points(successful_iterations, n_sample_points)
  } else {
    sample_points <- c(1, successful_iterations)
  }
  
  # standardize starting energy to "1" 
  energy_tmp <- res_min_conf$energy[[2]][sample_points]
  energy_before <- energy_tmp 
  if (res_min_conf$energy[[2]][1] > 0){
    energy_tmp <- res_min_conf$energy[[2]][sample_points] / res_min_conf$energy[[2]][1]
  }
  
  iteration_tmp <- res_min_conf$energy[[1]][sample_points]
  
  result <- tibble::tibble(sites = n_sites,
                           gamma = total_gamma,
                           mean_alpha = mean_alpha,
                           species = richness_vector,
                           richness_sd = richness_sd,
                           max_runs = max_runs,
                           energy_threshold = energy_threshold, 
                           energy_scaled = energy_tmp,
                           iteration = iteration_tmp,
                           energy_before = energy_before,
                           min_conf_t = time_min_conf,
                           #tabu = tabu,
                           replicate = replicate)
  
  
  # Write output
  if(isTRUE(writeRDS)) {
    saveRDS(result, file.path(paste0("results/spectre_", siminputrow, ".rds")))
  }
  
  return(result)
}

### 


load("vignettes/p.rda")

x <- foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
  # library(spectre)
  x <- sim_fun(REPLICATE,
                       p,
                       TRUE,
                       TRUE)
  return(1)
  
  
} 


