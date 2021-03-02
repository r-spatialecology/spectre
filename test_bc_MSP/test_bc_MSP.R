### test whether optimizing using bray curtis (bc) dissimilarity or commonness will be better. 
### Sebastian`s idea was that bc contains more biologically relevant information 

### use the Bangalore data... 

# Evaluate the spectre workflow using 36 farm sites and Bangalore data
# we will explore how much the algorithm, and each of the diversity models contributes to the overall error

library("spectre")
library("foreach")

# set parameters 
replicate <- 1:15 # replicates per parameter combination
max_runs <-100000 # 50000 
energy_threshold <- 0.0
verbose = TRUE

p <- tidyr::crossing(replicate, 
                     max_runs,
                     energy_threshold, 
                     verbose
)

p$siminputrow <- 1:dim(p)[1]


# Read in observed 2016 Bangalore bird species list 
# =================================================
data_list <- readRDS("./test_bc_MSP/BEN_2016_bird_obs.rds") 
data <- as.matrix(data_list)
data <- t(data) # transform dataset 
print(paste0(dim(data)[1], " species, ",  dim(data)[2], " sites")) # just a check

# Generate alpha_list, total_gamma &  "real" commonness matrix from data
n_sites <- dim(data)[2]
n_species <- dim(data)[1]
(total_gamma <- sum(rowSums(data) > 0))
observed_commonness <- spectre:::calculate_solution_commonness_rcpp(data)
mean_observed_commonness <- mean(abs(observed_commonness), na.rm = TRUE)
 # observed_sorensen <- spectre:::calculate_solution_sorensen_rcpp(data)

# Alpha diversity 
# ===============
(alpha_list <- colSums(data)) # original data

# from lm (full model)
alpha_model_full <- readRDS("./test_bc_MSP/BEN_richness_2016.rds")

# Beta-diversity / commonness
# ===========================

gdm_predictions <- readRDS("./test_bc_MSP/BEN_beta_2016.rds") # predicted Sorensen estimates, from Gabriel
gdm_predictions <- unlist(gdm_predictions) 

predicted_commonness <- spectre:::generate_commonness_matrix_from_gdm(gdm_predictions = gdm_predictions, 
                                                                      alpha_list = alpha_model_full)

# how far does the predicted commonness matrix deviate from observed data?
(MAE_target_observed <- mean(abs(predicted_commonness - observed_commonness), na.rm = TRUE))
(RCE_target_observed <- MAE_target_observed / mean_observed_commonness * 100)

### Predict alternative target matrices from a) predicted alpha and observed beta diversity
### and b) observed alpha and predicted beta-diversity. Both targets are evaluated against the observed data  

gdm_sorensen <- spectre:::generate_sorensen_matrix_from_gdm(gdm_predictions = gdm_predictions, 
                                                            n_sites = n_sites)


sim_fun <- function(siminputrow, parameters, writeRDS, verbose)  
{
  # Read and set parameters 
  p <- parameters[siminputrow, ]
  if(!( dim(p)[1] == 1)) {
    p <- p[1, ] # for quick testing 
  }
  energy_threshold <- p$energy_threshold
  replicate <- p$replicate
  max_runs <- p$max_runs
  # n_sample_points <- p$n_sample_points
  seed <- p$siminputrow
  
  # Set seed
  set.seed(seed)
  
  res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_model_full, # alpha_list, alpha_model_full, alpha_model_poly
                                                       total_gamma = total_gamma, 
                                                       target = predicted_commonness, # observed_commonness, predicted_commonness
                                                       max_iterations = max_runs,
                                                       # patience = 2500, 
                                                       # energy_threshold = energy_threshold,
                                                       verbose = verbose,
                                                       seed = seed) # 
  
  
  solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid)
  
  # save error between solution and (1) predicted commonness and (2) observed commonness in tibble
  # same for sorensen... 
  
  (MAE_sol_obs <- mean(abs(solution_commonness - observed_commonness), na.rm = TRUE)  )
  RCE_sol_obs <- MAE_sol_obs / mean_observed_commonness * 100
  
  (MAE_sol_pred <- mean(abs(solution_commonness - predicted_commonness), na.rm = TRUE)  )
  RCE_sol_pred <- MAE_sol_pred / mean_observed_commonness * 100
  
  
  r_1 <- tibble::tibble(replicate = replicate,
                        MAE_target_observed = MAE_target_observed, 
                        RCE_target_observed = RCE_target_observed,
                        MAE_sol_obs = MAE_sol_obs,
                        RCE_sol_obs = RCE_sol_obs,
                        MAE_sol_pred = MAE_sol_pred, 
                        RCE_sol_pred = RCE_sol_pred)
  saveRDS(r_1, file = paste0("./test_bc_MSP/res_bc/comm", siminputrow, ".rds")) # change this later! 
}

foreach(REPLICATE = 1:dim(p)[1], .export = c("p"), .packages = c("spectre")) %do% {
  print(paste0("Replicate: ", REPLICATE))
  
  sim_fun(REPLICATE,
          p,
          TRUE,
          TRUE)
  return(1)
} 




