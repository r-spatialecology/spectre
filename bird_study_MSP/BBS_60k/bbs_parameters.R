# Generate parameter table... 
# optimized for min_conf 

setwd("~/spectre/bird_study_MSP/BBS_60k")

replicate <- 1:3 # replicates per parameter combination
n_sites <- c(30, 40, 49) # number of sites  c(20, 40, 60)
n_species <- c(30, 50, 100, 150) # 30, 50, 100, 150) # 25, 50, 100, 150
total_gamma <- c(15) # evaluated gamma diversities  c(30, 60, 90) 
richness_vector <- c(0.1) # evaluated fractions of gamma diversity (mean richness) c(0.2, 0.4, 0.6) 
richness_sd <- 0 
max_runs <- 60000
energy_threshold <- 0.0
n_sample_points <- 50
tabu_percent <- c(0) # percent! 
verbose = TRUE

p <- tidyr::crossing(replicate, 
                     n_sites,
                     n_species,
                     total_gamma,
                     richness_vector,
                     richness_sd,
                     max_runs,
                     energy_threshold, 
                     n_sample_points,
                     tabu_percent,
                     verbose
)

p$mean_alpha <- round(p$total_gamma * p$richness_vector)
p$siminputrow <- 1:dim(p)[1]
save(p, file = "p.rda")
print(paste0(dim(p)[1], " parameter combinations will be tested"))

