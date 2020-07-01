# Generate parameter table... 
# optimized for min_conf 
replicate <- 1:25 # replicates per parameter combination
n_sites <- c(20, 40, 60) # number of sites  c(20, 40, 60)
total_gamma <- c(30, 60, 90) # evaluated gamma diversities  c(30, 60, 90) 
richness_vector <- c(0.1, 0.2, 0.4) # evaluated fractions of gamma diversity (mean richness) c(0.2, 0.4, 0.6) 
richness_sd <- 0 
max_runs <- 15000
energy_threshold <- 0.0
n_sample_points <- 50 
tabu_percent <- c(0) # percent! 

p <- tidyr::crossing(replicate, 
                     n_sites,
                     total_gamma,
                     richness_vector,
                     richness_sd,
                     max_runs,
                     energy_threshold, 
                     n_sample_points,
                     tabu_percent
)

p$mean_alpha <- round(p$total_gamma * p$richness_vector)
p$siminputrow <- 1:dim(p)[1]
save(p, file = "vignettes/p.rda")
print(paste0(dim(p)[1], " parameter combinations will be tested"))

      