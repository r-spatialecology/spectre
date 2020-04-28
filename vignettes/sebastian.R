library(spectre)

# welcome to the section of simulated (sim) data 
# ==============================================
# the aim is to test/compare algorithm performance under different conditions:
# i.e. number of sites, gamma diversity and average richness per site
# as a start, we assume no correlations between sites
total_gamma_sim <- 70 
n_sites_sim <- 25 

# average number of species per site, if equal richness across all sites is demanded, change "sd_sim" to "0".
# to allow both for constant and varying richness across sites, I used round( rnorm() ) instead of rpois() here. 
mean_alpha_sim <- 30 
sd_sim <- 10 # variation in mean / change to 0 if no variation is needed 
alpha_list_sim <- round(rnorm(n = n_sites_sim, mean = mean_alpha_sim, sd = sd_sim))
alpha_list_sim[alpha_list_sim < 1] <- 1 # only positive species numbers allowed 
print(alpha_list_sim) # check
hist(alpha_list_sim) # visual check of richness per site

# create simulated target
target_matrix_sim <- spectre:::generate_data_simple(total_gamma = total_gamma_sim, n_sites = n_sites_sim, alpha_list = alpha_list_sim)

# solve / optimize 
res_sim <- spectre:::optimizer(alpha_list_sim, total_gamma_sim, target_matrix_sim, 20000 )
spectre::plot_energy(res_sim)
spectre::plot_commonness(res_sim, target_matrix_sim)

# end of simulated data section 
# =============================





load("data/alpha_list_test.rda")
load("data/target_matrix_test.rda")
load("data/total_gamma_test.rda")

res1 <- spectre:::optimizer(alpha_list_test, total_gamma_test, target_matrix_test, 2000)

spectre::plot_energy(res1)
spectre::plot_commonness(res1, target_matrix_test)
res1$optimized_grid


load("data/alpha_list.rda")
load("data/target_matrix.rda")
load("data/estimated_gamma.rda")
# 15:30 - 16:10
res2 <- spectre:::optimizer(alpha_list, estimated_gamma, target_matrix, 100000) # ~20%

spectre::plot_energy(res2)
spectre::plot_commonness(res2, target_matrix)

res_spec <- c()
SITES <- 100
load("data/random_solution.rda")
for (SPECIES in 20:139) {
  
  current_solution_15 <- current_solution[1:SPECIES,1:SITES]
  alpha_15 <- vector()
  for (i in 1:SITES) {
    alpha_15[i] <- sum(current_solution_15[,i])
  }
  target_fake <- spectre:::calculate_solution_commonness_rcpp(current_solution_15)
 # saveRDS(alpha_15, file = "data/alpha100x20.rds")
 #  saveRDS(target_fake, file = "data/target100x20.rds")
  #res_15 <- spectre:::optimizer(alpha_15, SPECIES, target_fake, 0) # Optimization finished after 1633000180 steps w 35 species
  res_spec <- c(res_spec, spectre:::optimizer(alpha_15, SPECIES, target_fake, 1000000000)) # successful with 80 species
}
energy <- c()
for (i in seq(2,240,2)) {
  energy <- c(energy,min(res_spec[[i]]))
}

energy <- tibble(n_species = c(20:139),
                 energy = energy)

ggplot(energy, aes(x = n_species, y = energy)) +
  geom_point()

spectre::plot_commonness(res_15, target_fake)
spectre::plot_commonness(res_mc, target_fake)
grid15 <- res_15$optimized_grid

big_grid <- cbind(res2$optimized_grid, res2$optimized_grid, res2$optimized_grid, res2$optimized_grid, res2$optimized_grid)
big_alpha <- cbind(alpha_list,alpha_list,alpha_list,alpha_list,alpha_list)

target_fake <- spectre:::calculate_solution_commonness_rcpp(big_grid)
res2_fake <- spectre:::optimizer(big_alpha, estimated_gamma, target_fake, 600) # ~13%

target_fake <- spectre:::calculate_solution_commonness_rcpp(res2$optimized_grid)
res2_fake <- spectre:::optimizer(alpha_list, estimated_gamma, target_fake, 500)
spectre::plot_energy(res2_fake)
spectre::plot_commonness(res2_fake, target_fake)

res_list <- list()

seq <- seq(2,500,10)


for (i in seq(2,500,10)) {
  # smaller
  smaller <- 1:i
  alpha_list_small <- big_alpha[smaller]
  target_small <- target_fake[smaller, smaller]
  res_list[[i]] <- spectre:::optimizer(alpha_list_small, estimated_gamma, target_small, 600)
}

saveRDS(res_list, file = "results/res_big_list.rds")
saveRDS(energy_data, file = "results/energy_fake_sorted.rds")

energy_data <- tibble(n_sites = seq(2,500,10),
                      energy = seq(2,500,10))

for (i in 1:length(seq)) {
  energy_data$energy[i] <- min(res_list[[seq[i]-1]]$energy)
}

energy_data <- rbind(energy_data, c(500, 0.12727))

ggplot(energy_data, aes(x = n_sites, y = energy)) +
  geom_point()

spectre::plot_energy(res_fake_small)
spectre::plot_commonness(res_fake_small, target_fake_small)

res_old <- spectre::run_optimization(alpha_list, estimated_gamma, target_matrix, annealing = 0.01, max_runs = 200000, energy_threshold = 0.0, patience = 1000)
# ~40%

spectre::plot_energy(res_old)
spectre::plot_commonness(res_old, target_matrix)

# ~ 1.8 mins @ufom-p23
bench::mark(
  res2 <- spectre:::optimizer(alpha_list, estimated_gamma, target_matrix, 1000)
)
