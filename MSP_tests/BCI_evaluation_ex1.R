# Illustrate performance of the spectre algorithm/package using subsets of the Barro Colorado tree data 
# Condit R., Perez, R., Aguilar, S., Lao, S., Foster, R., Hubbell, S.P. 2019. 
# Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years, 2019 version.   
# https://doi.org/10.15146/5xcp-0d46.

# Create error ~ iteration plot data 

out_dir <- "MSP_tests/error_iteration"
getwd()
dir.create(file.path(out_dir))

library("spectre")
library("foreach")

# set parameters 
replicate <- 1 # replicates per parameter combination
n_sites <- c(8) # how large is the BCI subset used
n_species <- c(10) # how many species are sampled from subset of the BCI data 
max_iterations <- c(20) # how many iterations of the algorithm
energy_threshold <- 0.0
verbose = TRUE

parameters <- tidyr::crossing(replicate, 
                              n_sites,
                              n_species,
                              max_iterations,
                              energy_threshold, 
                              verbose)


parameters$siminputrow <- 5 # 1:dim(parameters)[1]

# load BCI data
data <- readRDS("data/BCI_tree8_MSP.rds") 
print(paste0(dim(data)[1], " species, ",  dim(data)[2], " sites"))


p <- parameters[1, ] # only used for quick testing 

(n_sites <- p$n_sites)
(n_species <- p$n_species) 
replicate <- p$replicate
max_iterations <- p$max_iterations
energy_threshold <- p$energy_threshold 
verbose <- p$verbose
siminputrow <- p$siminputrow
seed <- siminputrow
set.seed(seed)

# create BCI subsample
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

# check number of species in subsample 
available_species <- which(rowSums(temp_data) > 0)
n_available_species <- length(available_species)

# if (n_available_species >= n_species){

# check that sampled sites have a commonness > 1 (very small samples could have only sites that share no species)
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

# get crucial information from sampled data
(alpha_list <- colSums(sampled_data))
(total_gamma <- sum(rowSums(sampled_data) > 0))
time_before <- Sys.time()
res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_list, 
                                                   total_gamma = total_gamma, 
                                                   target = target_commonness, 
                                                   max_iterations = max_iterations,
                                                   # patience = 2500, 
                                                   #energy_threshold = energy_threshold,
                                                   verbose = verbose,
                                                   seed = seed) #



### Evaluate error between target and solution
(solution_commonness <- spectre:::calculate_solution_commonness_rcpp(res_min_conf$optimized_grid))

target_commonness

# Absolute error, calculated from solution grid
sum(abs(target_commonness - solution_commonness), na.rm = TRUE)

# "Absolute error" from spectre output
res_min_conf$error$error[length(res_min_conf$error$error)] # after last iteration
res_min_conf$error$error # all iterations 


# I would expect that error from optimized solution and spectre output error should match... 

