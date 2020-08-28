# Test min_conf_algorithms with BCI data

setwd("~/spectre/BCI_tidy")

library("spectre")
library("foreach")

load("p_BCI.rda")
p

data <- readRDS("data/BCI_tree8_MSP.rds") 

dim(data) 
print(paste0(dim(data)[1], " species, ",  dim(data)[2], " sites"))

sim_fun <- function(siminputrow, parameters, writeRDS, verbose)  # Code mostly stolen from Jan Salecker... 
{
  # Read and set parameters 
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) {
    p <- p[1, ]
  }
  (n_sites <- p$n_sites)
  (n_species <- p$n_species) 
  replicate <- p$replicate
  max_runs <- p$max_runs
  n_sample_points <- p$n_sample_points
  energy_threshold <- p$energy_threshold 
  verbose <- p$verbose
  siminputrow <- p$siminputrow
  seed <- siminputrow
  
  set.seed(seed)
  
  # create data subsample only with sampled sites 
  
  sampled_data <- matrix(nrow = n_species, ncol = n_sites, data = 0)
  
  if(n_sites > dim(data)[2]) { # sites in columns 
    n_sites <- dim(data)[2] 
    print("Could not sample all requested sites: you reached maximum number of sites... ")
  }
  sampled_sites <-  sample(dim(data)[2], n_sites, replace = FALSE) 
  
  if (n_species > dim(data)[1]) { # species in rows
    n_species <- dim(data)[1]
    print("Could not sample all requested species: You reached maximum number of species ... ")
  }
  
  temp_data <- data[, sampled_sites]
  
  # check whether there are enough species in subsample 
  
  available_species <- which(rowSums(temp_data) > 0)
  n_available_species <- length(available_species)
  if (n_available_species >= n_species){
    
    
    sum_commonness <- 0
    while (sum_commonness < 1) {
      sampled_species <-  sample(available_species, n_species, replace = FALSE) 
      sampled_data[1:n_species, ] <- temp_data[sampled_species, ]
      
      # check commonness of sampled sites
      (target_commonness <- spectre:::calculate_solution_commonness_rcpp( sampled_data ) )
      sum_commonness <- sum(target_commonness, na.rm = TRUE)
      print(paste0("Sum commonness is = ", sum_commonness))
    }
    
    print(paste0(dim(sampled_data)[1], " species, ",  dim(sampled_data)[2], " sites"))
    
    # get & save crucial information from sampled data
    
    (alpha_list <- colSums(sampled_data))
    (total_gamma <- sum(rowSums(sampled_data) > 0))
    mean_richness <- mean(alpha_list)
    mean_commonness <- mean(abs(target_commonness), na.rm = TRUE)
    
    start_time <- Sys.time()
    res_min_conf <- spectre::run_optimization_min_conf_0(alpha_list = alpha_list, 
                                                         total_gamma = total_gamma, 
                                                         target = target_commonness, 
                                                         max_iterations = max_runs,
                                                         patience = 2500, 
                                                         energy_threshold = energy_threshold,
                                                         verbose = verbose,
                                                         seed = seed) # 
    
    end_time <- Sys.time()
    time_min_conf <- as.numeric(end_time - start_time)
    
    solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
    
    x <- table(abs(solution_commonness - target_commonness))
    
    mean_real_commonness_error <- mean(abs(solution_commonness - target_commonness), na.rm = TRUE) 
    # plot(x)
    x_df <- data.frame(x)
    names(x_df) <- c("Commonness", "Count")
    
    r_1 <- tibble::tibble(x_df, 
                          replicate = replicate,
                          mean_commonness = mean_commonness,
                          mean_real_commonness_error = mean_real_commonness_error,
                          final_distance_D = res_min_conf$energy[[2]][length(res_min_conf$energy[[2]])],
                          n_sites = n_sites,
                          n_species = n_species,
                          mean_richness = mean_richness,
                          iterations = max_runs,
                          min_conf_t = time_min_conf)
    saveRDS(r_1, file = paste0("./table_res/comm", siminputrow, ".rds")) 
    
    
    
    
    
    # only keep a fraction of results (first, sample_points in between , last)
    successful_iterations <- sum(!is.na(res_min_conf$energy[2]))
    
    if(successful_iterations >= n_sample_points){
      sample_points <- spectre:::get_sample_points(successful_iterations, n_sample_points)
    } else {
      sample_points <- c(1, successful_iterations)
    }
    
    iteration_tmp <- res_min_conf$energy[[1]][sample_points]
    energy_before <- res_min_conf$energy[[2]][sample_points]
    
    result <- tibble::tibble(sites = n_sites,
                             species = n_species,
                             max_runs = max_runs,
                             energy_threshold = energy_threshold, 
                             mean_real_commonness_error = mean_real_commonness_error,  
                             iteration = iteration_tmp,
                             energy_before = energy_before,
                             min_conf_t = time_min_conf,
                             replicate = replicate)
    
    
    # Write output
    if(isTRUE(writeRDS)) {
      saveRDS(result, file.path(paste0("res/", siminputrow, ".rds"))) 
    }
    
    return(1)
  }
}
### 



foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
  print(paste0("Replicate: ", REPLICATE))
  
  sim_fun(REPLICATE,
          p,
          TRUE,
          TRUE)
  return(1)
} 




