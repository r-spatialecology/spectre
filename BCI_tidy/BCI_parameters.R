# Generate parameter table... 
# optimized for min_conf 

setwd("~/spectre/BCI_tidy")

replicate <- 1 # replicates per parameter combination
n_sites <- c(20, 40, 60, 80, 100) # number of sites  c(20, 40, 60)
n_species <- c(30, 50 , 100, 150) # 25, 50, 100, 150 50, 100, 150
max_runs <- c(2000, 10000, 20000)
energy_threshold <- 0.0
n_sample_points <- 500
verbose = TRUE

p <- tidyr::crossing(replicate, 
                     n_sites,
                     n_species,
                     max_runs,
                     energy_threshold, 
                     n_sample_points,
                     verbose
)


p$siminputrow <- 1:dim(p)[1]

save(p, file = "p_BCI.rda")
print(paste0(dim(p)[1], " parameter combinations will be tested"))

