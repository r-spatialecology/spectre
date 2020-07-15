# Test min_conf_algorithms with bird breeding survey (BBS) data

setwd("~/spectre/bird_study_MSP/BBS_60k")

library("spectre")
library("foreach")
library("betapart")

load("p.rda")


data(bbsData)
data <- bbs1980

dim(data)
data <- t(data) # spectre expects species in rows, sites in columns 
print(paste0(dim(data)[1], " species, ",  dim(data)[2], " sites"))

sim_fun <- function(siminputrow, parameters, writeRDS)  # Code mostly stolen from Jan Salecker... 
{
  # Read and set parameters 
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) {
    p <- p[2, ]
  }
  (n_sites <- p$n_sites)
  (n_species <- p$n_species) 
  mean_alpha <- p$mean_alpha 
  richness_vector <- p$richness_vector
  richness_sd <- p$richness_sd 
  replicate <- p$replicate
  max_runs <- p$max_runs
  seed <- p$siminputrow
  n_sample_points <- p$n_sample_points
  energy_threshold <- p$energy_threshold 
  tabu <- round(p$tabu_percent * p$n_sites / 100)
  verbose <- p$verbose  
  
  # Set seed
  set.seed(seed)
  
  # subsample data, check that subsample has sites with common species at all, 
  # bc spectre otherwise returns "nan" and the script crashes...
  
  sum_commonness <- 0 
  while(sum_commonness < 1){
    
    if (n_species > dim(data)[1]) { # species in rows
      n_species <- dim(data)[1]
      print("Could not sample all requested species: You reached maximum number of species ... ")
    }
    sampled_species <-  sample(dim(data)[1], n_species, replace = FALSE) 
    
    if(n_sites > dim(data)[2]) { # sites in columns 
      n_sites <- dim(data)[2] 
      print("Could not sample all requested sites: you reached maximum number of sites... ")
    }
    sampled_sites <-  sample(dim(data)[2], n_sites, replace = FALSE) 
    
    
    
    # target species x site list 
    species_list <- data[sampled_species, sampled_sites]
    dim(species_list)
    (alpha_list <- colSums(species_list))
    length(alpha_list) == n_sites
    (total_gamma <- sum(rowSums(species_list) > 0))
    
    
    
    if (total_gamma < n_species) {
      dropped_species <- which(rowSums(species_list) == 0)
      species_list <- species_list[- dropped_species, ]
      print(paste0("Realized gamma is: ",total_gamma, " n species is: ", n_species ) )
      print(paste0("Dimensions of species list: ", dim(species_list)[1]))
    }
    
    
    (target_commonness <- spectre:::calculate_solution_commonness_rcpp( species_list ) )# transpose !!!
    target_commonness[upper.tri(target_commonness, diag=TRUE)] <- NA
    class(target_commonness)
    sum_commonness <- sum(target_commonness, na.rm = TRUE)
    print(paste0("Sum commonness is = ", sum_commonness))
  }
  
  print(paste0("How many species sampled? A:", dim(species_list)[1]))
  
  #print("start")
  # Executing original algorithm:
  
  start_time <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf_0(alpha_list = alpha_list, 
                                                       total_gamma = total_gamma, 
                                                       target = target_commonness, 
                                                       max_iterations = max_runs,
                                                       patience = 2500, 
                                                       energy_threshold = energy_threshold,
                                                       verbose = verbose,
                                                       seed = seed) # 
  
  solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
  
  (correct_pairs <- sum((solution_commonness - target_commonness) == 0, na.rm = TRUE) )
  (incorrect_pairs <- sum((solution_commonness - target_commonness) != 0, na.rm = TRUE) )
  (correctly_predicted_site_pairs <- correct_pairs / (correct_pairs + incorrect_pairs))
  
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
  
  best_energy <- res_min_conf$best_energy # problematic! 
  
  result <- tibble::tibble(sites = n_sites,
                           species = n_species,
                           total_gamma = total_gamma, 
                           mean_alpha = mean_alpha,
                           # species = richness_vector,
                           richness_sd = richness_sd,
                           max_runs = max_runs,
                           energy_threshold = energy_threshold, 
                           energy_scaled = energy_tmp,
                           iteration = iteration_tmp,
                           energy_before = energy_before,
                           min_conf_t = time_min_conf,
                           best_energy, 
                           correct_pairs = correctly_predicted_site_pairs, 
                           replicate = replicate)
  
  
  # Write output
  if(isTRUE(writeRDS)) {
    saveRDS(result, file.path(paste0(siminputrow, ".rds"))) 
  }
  
  return(1)
}

### 



foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
  print(paste0("Replicate: ", REPLICATE))
  
  sim_fun(REPLICATE,
          p,
          TRUE)
  return(1)
} 




