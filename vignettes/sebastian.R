library(spectre)

# welcome to the section of simulated (sim) data 
# ==============================================
# the aim is to test/compare algorithm performance under different conditions:
# i.e. number of sites, gamma diversity and average richness per site
# as a start, we assume no correlations between sites
total_gamma_sim <- 30
n_sites_sim <- 25 

# average number of species per site, if equal richness across all sites is demanded, change "sd_sim" to "0".
# to allow both for constant and varying richness across sites, I used round( rnorm() ) instead of rpois() here. 
mean_alpha_sim <- 12 
sd_sim <- 5 # variation in mean / change to 0 if no variation is needed 
alpha_list_sim <- round(rnorm(n = n_sites_sim, mean = mean_alpha_sim, sd = sd_sim))
alpha_list_sim[alpha_list_sim < 1] <- 1 # only positive species numbers allowed 
print(alpha_list_sim) # check
hist(alpha_list_sim) # visual check of richness per site

# create simulated target
target_matrix_sim <- spectre:::generate_data_simple(total_gamma = total_gamma_sim, n_sites = n_sites_sim, alpha_list = alpha_list_sim)
fixed_species <- matrix()

measure <- tibble::tibble(
  measure = vector(length = 400),
  n =  as.factor(c(rep(1, 100), rep(10, 100), rep(100, 100), rep(1000, 100))))
for (i in 1:100) {
  measure$measure[i] <- spectre:::calc_random_energy(1, alpha_list_sim, total_gamma_sim, target_matrix_sim, i, "max")
  measure$measure[i + 100] <- spectre:::calc_random_energy(10, alpha_list_sim, total_gamma_sim, target_matrix_sim, i, "max")
  measure$measure[i + 200] <- spectre:::calc_random_energy(100, alpha_list_sim, total_gamma_sim, target_matrix_sim, i, "max")
  measure$measure[i + 300] <- spectre:::calc_random_energy(1000, alpha_list_sim, total_gamma_sim, target_matrix_sim, i, "max")
}

library(ggplot2)

p <- ggplot(measure, aes(x = n, y = measure))
p + geom_boxplot() + 
  labs(title = "Averaged measure for a random solution (max norm)", x = "# repetitions")
ggsave("random_max.png")


res_sim0_sum <- run_optimization_min_conf_0(alpha_list = alpha_list_sim,
                                            total_gamma = total_gamma_sim,
                                            target = target_matrix_sim,
                                            fixed_species = NULL,
                                            partial_solution = NULL,
                                            max_iterations = 15000,
                                            energy_threshold = 0.0,
                                            seed = 1,
                                            verbose = TRUE,
                                            norm = "sum")

res_sim0_euclid <- run_optimization_min_conf_0(alpha_list = alpha_list_sim,
                                               total_gamma = total_gamma_sim,
                                               target = target_matrix_sim,
                                               fixed_species = NULL,
                                               partial_solution = NULL,
                                               max_iterations = 15000,
                                               energy_threshold = 0.0,
                                               seed = 1,
                                               verbose = TRUE,
                                               norm = "euclid")

res_sim0_max <- run_optimization_min_conf_0(alpha_list = alpha_list_sim,
                                            total_gamma = total_gamma_sim,
                                            target = target_matrix_sim,
                                            fixed_species = NULL,
                                            partial_solution = NULL,
                                            max_iterations = 15000,
                                            energy_threshold = 0.0,
                                            seed = 1,
                                            verbose = TRUE,
                                            norm = "max")
plot_energy(res_sim0_sum)
spectre:::calc_energy(spectre:::calculate_solution_commonness_rcpp(res_sim0_sum$optimized_grid), target_matrix_sim)
plot_commonness(res_sim0_sum, target_matrix_sim)
plot_energy(res_sim0_euclid)
spectre:::calc_energy(spectre:::calculate_solution_commonness_rcpp(res_sim0_euclid$optimized_grid), target_matrix_sim)
plot_commonness(res_sim0_euclid, target_matrix_sim)
plot_energy(res_sim0_max)
spectre:::calc_energy(spectre:::calculate_solution_commonness_rcpp(res_sim0_max$optimized_grid), target_matrix_sim)
plot_commonness(res_sim0_max, target_matrix_sim)




true_solution <- matrix(data = 0, 
                        nrow = total_gamma_sim, 
                        ncol = n_sites_sim)

for (n in seq_len(n_sites_sim)) {
  
  change_locations <- spectre::rcpp_sample(x = seq_len(total_gamma_sim),
                                           size = alpha_list_sim[n])
  
  true_solution[change_locations, n] <- 1
}

target_matrix_sim <- spectre:::calculate_solution_commonness_rcpp(true_solution)

fixed_species <- matrix(data = c(c(0,1,1,1,0,0), rep(0, 54)), nrow = 6, ncol = 10)
fixed_species <- true_solution


# solve / optimize 
energy <- 0
energy0 <- 0
energy1 <- 0
energy2 <- 0
for (i in 1:10) {
  res_sim <- spectre:::optimizer_min_conf(alpha_list_sim, 
                                          total_gamma_sim, 
                                          target_matrix_sim, 5000, 0.0)
  energy <- energy + min(res_sim$energy$energy)
  
  res_sim0 <- spectre:::optimizer_min_conf0(alpha_list = alpha_list_sim, 
                                            total_gamma = total_gamma_sim, 
                                            target = target_matrix_sim,
                                            fixed_species = fixed_species,
                                            max_iterations = 5000, 
                                            energy_threshold = 0.0)
  energy0 <- energy0 + min(res_sim0$energy$energy)
  
  
  
  res_sim1 <- spectre:::optimizer_min_conf1(alpha_list = alpha_list_sim, 
                                            total_gamma = total_gamma_sim, 
                                            target = target_matrix_sim,
                                            fixed_species = fixed_species,
                                            tabu = 0,
                                            fixed_partial_solution = TRUE,
                                            max_iterations = 5000, 
                                            energy_threshold = 0.0)
  energy1 <- energy1 + min(res_sim1$energy$energy)
  
  res_sim2 <- spectre:::optimizer_min_conf2(alpha_list = alpha_list_sim, 
                                            total_gamma = total_gamma_sim, 
                                            target = target_matrix_sim,
                                            fixed_species = fixed_species,
                                            tabu = 0,
                                            fixed_partial_solution = TRUE,
                                            max_iterations = 5000, 
                                            energy_threshold = 0.0)
  energy2 <- energy2 + min(res_sim2$energy$energy)
}

energy0 <- 0
energy1 <- 0
for (i in 1:10) {
  res_sim0 <- spectre:::optimizer_min_conf0(alpha_list = alpha_list_sim, 
                                            total_gamma = total_gamma_sim, 
                                            target = target_matrix_sim,
                                            fixed_species = fixed_species,
                                            tabu = 0,
                                            fixed_partial_solution = TRUE,
                                            max_iterations = 25000, 
                                            energy_threshold = 0.0)
  energy0 <- energy0 + min(res_sim0$energy$energy) / 10
  
  res_sim1 <- spectre:::optimizer_min_conf1(alpha_list = alpha_list_sim, 
                                            total_gamma = total_gamma_sim, 
                                            target = target_matrix_sim,
                                            fixed_species = fixed_species,
                                            tabu = 0,
                                            fixed_partial_solution = TRUE,
                                            max_iterations = 25000, 
                                            energy_threshold = 0.0)
  energy1 <- energy1 + min(res_sim1$energy$energy) / 10
}

saveRDS(res_sim1, file = "./res1.RDS")

res_sim2$optimized_grid

spectre::plot_energy(res_sim)
spectre::plot_energy(res_sim0)
spectre::plot_energy(res_sim01)
spectre::plot_energy(res_sim2)

spectre::plot_commonness(res_sim0, target_matrix_sim)
spectre::plot_commonness(res_sim01, target_matrix_sim)
spectre::plot_commonness(res_sim2, target_matrix_sim)
spectre::plot_commonness(res_sim3, target_matrix_sim)

alpha_list_res1 <- alpha_list_sim
alpha_list_res2 <- alpha_list_sim
for (site in 1:n_sites_sim) {
  alpha_list_res1[site] <-  sum(res_sim1$optimized_grid[,site])
  alpha_list_res2[site] <-  sum(res_sim2$optimized_grid[,site])
}
res_sim2$optimized_grid[,1]

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
