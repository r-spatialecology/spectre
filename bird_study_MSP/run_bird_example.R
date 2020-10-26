# Test min_conf_algorithms 
###todo CHECK SAMPLE POINTS FUNCTION!!! 
library("spectre")
library("foreach")
library("betapart")

load("bird_study_MSP/p.rda")
p
data(bbsData)
data <- bbs1980

# Error in { : 
# task 7 failed - "cannot take a sample larger than the population when 'replace = FALSE'"

# Error in { : 
# task 1 failed - "Column name `species` must not be duplicated." from tibble::! 


### 
sim_fun <- function(siminputrow, parameters, writeRDS, verbose)  # Code mostly stolen from Jan Salecker... 
{
  # Read and set parameters p <- p[1, ]
  p <- parameters[siminputrow, ]
 # if(!( dim(p)[1] == 1)) {
  #  p <- p[1, ]
  #}
  n_sites <- p$n_sites
  n_species <- p$n_species 
  mean_alpha <- p$mean_alpha 
  richness_vector <- p$richness_vector
  richness_sd <- p$richness_sd 
  # total_gamma <- p$total_gamma
  replicate <- p$replicate
  max_runs <- p$max_runs
  n_sample_points <- p$n_sample_points
  energy_threshold <- p$energy_threshold 
  tabu <- round(p$tabu_percent * p$n_sites / 100)
  
  # Set random seed
  set.seed(replicate)
  
  nobirds <- FALSE
  if (nobirds){
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
    
  }
  
  ### subsample bird data
  sum_commonness <- 0
  while(sum_commonness < 1){
    if(n_sites > dim(data)[1]) {
      n_sites <- dim(data)[1] 
      print("Could not sample all requested sites: you reached maximum number of sites... ")
    }
    sampled_sites <-  sample(dim(data)[1], n_sites, replace = FALSE) 
    if (n_species > dim(data)[2]) {
      n_species <- dim(data)[2]
      print("Could not sample all requested species: You reached maximum number of species ... ")
    }
    sampled_species <-  sample(dim(data)[2], n_species, replace = FALSE) 
    species_list <- data[sampled_sites, sampled_species ]
    dim(species_list)
    (alpha_list <- rowSums(species_list))
    (total_gamma <- sum(colSums(species_list) > 0))
    # (gamma <- betapart::betapart.core(species_list)$St)
    total_gamma <- as.numeric(total_gamma)
    # species_list <- t(as.matrix(species_list))
    
    species_list <- t(species_list)
    
    (target_commonness <- spectre:::calculate_solution_commonness_rcpp( species_list ) )# transpose !!!
    target_commonness[upper.tri(target_commonness, diag=TRUE)] <- NA
    class(target_commonness)
    sum_commonness <- sum(target_commonness, na.rm = TRUE)
    print(paste0("Sum commonness is = ", sum_commonness))
  }
  
  #print("start")
  # Executing original algorithm:
  
  start_time <- Sys.time()
  res_min_conf <- spectre::run_optimization_min_conf_0(alpha_list = alpha_list, 
                                                       total_gamma = total_gamma, 
                                                       target = target_commonness, 
                                                       max_iterations = max_runs,
                                                       patience = 2500, 
                                                       energy_threshold = energy_threshold,
                                                       # tabu = tabu,
                                                       verbose = verbose) # 
  
  solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
  
  (correct_pairs <- sum((solution_commonness - target_commonness) == 0, na.rm = TRUE) )
  (incorrect_pairs <- sum((solution_commonness - target_commonness) != 0, na.rm = TRUE) )
  (correctly_predicted_site_pairs <- correct_pairs / (correct_pairs + incorrect_pairs))
  # print(paste0("Correctly predicted pairs: ", correctly_predicted_site_pairs * 100 , " % \n"))
  
  
  plot(res_min_conf$energy[[1]], res_min_conf$energy[[2]])
  abline(h=res_min_conf$best_energy)
  
  # calculate proportion of correctly predicted pairs with correct commonness
  
  
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
                           correct_pairs = correctly_predicted_site_pairs * 100, 
                           #tabu = tabu,
                           replicate = replicate)
  
  
  # Write output
  if(isTRUE(writeRDS)) {
    saveRDS(result, file.path(paste0("bird_study_MSP/results100/spectre_", siminputrow + 60, ".rds")))
  }
  
  return(result)
}

### 



# x <- foreach(REPLICATE = 1:2, .export = c("p"), .packages = c("spectre")) %do% {
x <- foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
  # library(spectre)
  print(paste0("Replicate: ", REPLICATE))
  x <- sim_fun(REPLICATE,
               p,
               TRUE,
               TRUE)
  return(1)
  
  
} 


