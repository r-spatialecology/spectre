# Generate parameter table... 
# optimized for min_conf 

setwd("~/spectre/BCI_tidy")

replicate <- 1:5 # 1:5 # replicates per parameter combination
n_sites <- c(30, 60, 90, 120, 150) # , 200, 250, 300) # c(20, 40, 60, 80, 100) # number of sites  c(20, 40, 60)
n_species <- c(25, 50, 100, 150, 200, 250, 300) # 
max_runs <- c(10000, 500000)
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

