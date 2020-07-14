# Test min_conf_algorithms 

setwd("~/spectre/bird_study_MSP/BCI_known_sites")

library("spectre")
library("foreach")
library("betapart")


load("p_BCI.rda")
p$siminputrow

data <- readRDS("~/spectre/data/BCI_tree8_MSP.rds")
dim(data)



### 
sim_fun <- function(siminputrow, parameters, writeRDS, verbose)  # Code mostly stolen from Jan Salecker... 
{
  # Read and set parameters p <- p[1, ]
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) { # only for easy testing, may be deleted later... 
    p <- p[1, ]
  }
  p
  n_sites <- p$n_sites
  known_sites <- p$known_sites
  n_species <- p$n_species 
  mean_alpha <- p$mean_alpha 
  richness_vector <- p$richness_vector
  richness_sd <- p$richness_sd 
  seed <- p$siminputrow
  max_runs <- p$max_runs
  n_sample_points <- p$n_sample_points
  energy_threshold <- p$energy_threshold 
  tabu <- round(p$tabu_percent * p$n_sites / 100)
  verbose <- p$verbose
  
  # Set random seed
  set.seed(seed)
  # subsample data, create species list, derive alpha diversity per site and 
  # gamma diversity from species list
  
  sum_commonness <- 0
  while(sum_commonness < 1){
    
    if(n_sites > dim(data)[2]) {
      n_sites <- dim(data)[2] 
      print("Could not sample all requested sites: you reached maximum number of sites... ")
    }
    sampled_sites <-  sample(dim(data)[2], n_sites, replace = FALSE) 
    
    if (n_species > dim(data)[1]) {
      n_species <- dim(data)[1]
      print("Could not sample all requested species: You reached maximum number of species ... ")
    }
    sampled_species <-  sample(dim(data)[1], n_species, replace = FALSE) 
    
    species_list <- data[sampled_species, sampled_sites]
    dim(species_list)
    (alpha_list <- colSums(species_list))
    (total_gamma <- sum(rowSums(species_list) > 0))
    
    # drop all species that do not occur in the sampled sites
    if (total_gamma < n_species) {
    dropped_species <- which(rowSums(species_list) == 0)
    species_list <- species_list[- dropped_species, ]
    }
    # total_gamma <- n_species # MSP: changed this line for the known species approach. 
    
    # total_gamma <- as.numeric(total_gamma)
    
    (target_commonness <- spectre:::calculate_solution_commonness_rcpp( species_list ) )# transpose !!!
    target_commonness[upper.tri(target_commonness, diag=TRUE)] <- NA
    class(target_commonness)
    sum_commonness <- sum(target_commonness, na.rm = TRUE)
    print(paste0("Sum commonness is = ", sum_commonness))
  }
  
  # chose fixed species
  fixed_species <- matrix(data = 0, nrow = total_gamma, ncol = n_sites)
  fixed_species[, 1:known_sites] <- species_list[, 1:known_sites]
  
  
  
  start_time <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf_0(alpha_list = alpha_list, 
                                                       total_gamma = total_gamma, 
                                                       target = target_commonness, 
                                                       fixed_species = fixed_species, 
                                                       max_iterations = max_runs,
                                                       patience = 2500, 
                                                       seed = seed,
                                                       energy_threshold = energy_threshold,
                                                       verbose = verbose) # 
  
  solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
  
  # correctly predicted commonness pairs
  
  (correct_pairs <- sum((solution_commonness - target_commonness) == 0, na.rm = TRUE) )
  (incorrect_pairs <- sum((solution_commonness - target_commonness) != 0, na.rm = TRUE) )
  (correctly_predicted_site_pairs <- correct_pairs / (correct_pairs + incorrect_pairs))
  # print(paste0("Correctly predicted pairs: ", correctly_predicted_site_pairs * 100 , " % \n"))
  
  # correctly predicted species 
  # like in Mokany et al. 2011 
  # number of observed species correctly predicted == observed in predicted 
  n_correctly_predicted <- sum(  which(res_min_conf$optimized_grid == 1) %in%  which(species_list == 1))
  
  n_observed_species <- length(which(species_list == 1))
  
  correctly_predicted_species <- n_correctly_predicted / n_observed_species 
  
  # print(paste0("Correctly predicted species: ", correctly_predicted_species ))
  
  # plot(res_min_conf$energy[[1]], res_min_conf$energy[[2]], ylim = c(0, 0.5), ylab = "energy", xlab = "iterations")
  
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
                           known_sites = known_sites, 
                           mean_alpha = mean_alpha,
                           # species = richness_vector,
                           richness_sd = richness_sd,
                           max_runs = max_runs,
                           energy_threshold = energy_threshold, 
                           energy_scaled = energy_tmp,
                           iteration = iteration_tmp,
                           energy_before = energy_before,
                           min_conf_t = time_min_conf,
                           # best_energy, 
                           correct_pairs = correctly_predicted_site_pairs * 100, 
                           correctly_predicted_species = correctly_predicted_species, 
                           #tabu = tabu,
                           replicate = seed)
  
  
  # Write output
  if(isTRUE(writeRDS)) {
    saveRDS(result, file.path(paste0(siminputrow, ".rds")))  
  }
  
  return(result)
}

### 




x <- foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
 
  print(paste0("Replicate: ", REPLICATE))
  x <- sim_fun(REPLICATE,
               p,
               TRUE,
               TRUE)
  return(1)
  
  
} 




