library(spectre)

# create parameter combinations:
# we evaluate all parameter combinations, so far without replicates
# mean richness per site is implemented as fraction of gamma diversity
# ======================================================================
n_sites <- c(25, 50) # number of sites  
gamma_vector <- c(40, 60, 80) # evaluated gamma diversities , 60, 80, 100, 120
richness_vector <- c(0.3, 0.5, 0.7) # evaluated fractions of gamma diversity (mean richness) , 0.5, 0.6

# parameter combinations are created here:
total_gamma <- rep(gamma_vector, each = length(richness_vector))
mean_alpha <- round(total_gamma * richness_vector) # average species per site 
sd_sim <- rep(0, length(mean_alpha)) # mean_alpha / 10 # variation in per-site richness / change to 0 if all sites in a run should have similar richness 

# we do not keep all calculated energy values, but the first, n_sample_points in between, and the last. 
# [note: some runs stop after 1 iteration (usually if richness is a high fraction of gamma diversity),
# thus a workaround to only keep one sampling point is found below...] 
n_sample_points <- 22 # how many sample points to keep 

# create lists to store results
RES_RND_SWITCH_ALL <- list() # "random switching"
RES_RND_SWITCH <- list()

RES_MIN_CONF_ALL <- list() # "8 queens"
RES_MIN_CONF <- list()

RES_BACKTRACKING_ALL <- list() # backtracking
RES_BACKTRACKING <- list()

N_SITES <- 1 # for quick testing
TOTAL_GAMMA <- 1

for (N_SITES in 1:length(n_sites)){
  for (TOTAL_GAMMA in 1:length(total_gamma)){
    
    alpha_list <- round(rnorm(n = n_sites[N_SITES], mean = mean_alpha[TOTAL_GAMMA], sd = sd_sim[TOTAL_GAMMA]))
    
    # check that richness per site is > 0 and <= gamma diversity. 
    # critical variables are "sd_sim" (variation in per site richness) and richness_vector 
    # (mean richness as fraction of gamma diversity)
    alpha_list[alpha_list < 1] <- 1 # only positive species numbers allowed 
    alpha_list[alpha_list > total_gamma[TOTAL_GAMMA]] <- total_gamma[TOTAL_GAMMA] #  
    
    # create target
    target_matrix_sim <- spectre:::generate_data_simple(total_gamma = total_gamma[TOTAL_GAMMA], n_sites = n_sites[N_SITES], alpha_list = alpha_list)
    
    # ====================================
    # solve / optimize using all algorithms
    # =====================================
    
    # random switching  
    # ================
    
    rd <- TRUE
    if(rd){
      max_run_rnd_switch <- 5000 # used also for result "sampling"
      res_rnd_switch <- spectre::run_optimization_rnd_switch(alpha_list = alpha_list, 
                                                             total_gamma = total_gamma[TOTAL_GAMMA], 
                                                             target = target_matrix_sim, 
                                                             max_runs = max_run_rnd_switch,
                                                             energy_threshold = 0.1,
                                                             patience = 500)
      
      # sample_points <- spectre:::get_sample_points(max_run_rnd_switch, n_sample_points)
      # # results from rnd_switch algorithm
      # sample_energy_rnd_switch <- res_rnd_switch$energy[[2]][sample_points]
      # sample_i_rnd_switch <- res_rnd_switch$energy[[1]][sample_points]
      
      ### copy section
      # only keep a fraction of results (first, sample_points in between , last)
      successful_iterations <- sum(!is.na(res_rnd_switch$energy[2]))
      print(paste0("Successful iterations of random_switching: ", successful_iterations))
      
      if(successful_iterations >= n_sample_points){
        sample_points <- spectre:::get_sample_points(successful_iterations, n_sample_points)
      } else {
        sample_points <- c(1, successful_iterations)
        print("Less successful iterations than sampling points. \n")
      }
      
      
      sites_tmp <- rep(n_sites[N_SITES], length(sample_points))
      total_gamma_tmp <- rep(total_gamma[TOTAL_GAMMA], length(sample_points))
      species_tmp <- rep(mean_alpha[TOTAL_GAMMA], length(sample_points))
      energy_tmp <- res_rnd_switch$energy[[2]][sample_points]
      iteration_tmp <- res_rnd_switch$energy[[1]][sample_points]
      
      RES_RND_SWITCH[[TOTAL_GAMMA]] <- cbind(sites_tmp,
                                             total_gamma_tmp, 
                                             species_tmp, 
                                             energy_tmp, 
                                             iteration_tmp)
    }
    
    # min-conflict (8-Queens)
    # =======================
    
    min_conf_iterations <- 5000 # used also for result "sampling" # 5000 default??? 
    res_min_conf <- spectre::run_optimization_min_conf(alpha_list = alpha_list, 
                                                       total_gamma = total_gamma[TOTAL_GAMMA], 
                                                       target = target_matrix_sim, 
                                                       max_runs = min_conf_iterations,
                                                       energy_threshold = 0.1) 
    
    # only keep a fraction of results (first, sample_points in between , last)
    successful_iterations <- sum(!is.na(res_min_conf$energy[2]))
    print(paste0("Successful iterations: ", successful_iterations))
    
    if(successful_iterations >= n_sample_points){
      sample_points <- spectre:::get_sample_points(successful_iterations, n_sample_points)
    } else {
      sample_points <- c(1, successful_iterations)
    }
    
    sites_tmp <- rep(n_sites[N_SITES], length(sample_points))
    total_gamma_tmp <- rep(total_gamma[TOTAL_GAMMA], length(sample_points))
    species_tmp <- rep(mean_alpha[TOTAL_GAMMA], length(sample_points))
    energy_tmp <- res_min_conf$energy[[2]][sample_points]
    iteration_tmp <- res_min_conf$energy[[1]][sample_points]
    
    RES_MIN_CONF[[TOTAL_GAMMA]] <- cbind(sites_tmp,
                                         total_gamma_tmp, 
                                         species_tmp, 
                                         energy_tmp, 
                                         iteration_tmp)
    
    # backtracking
    # ============
    # backtracking, in contrast to random switching and min-conf, gives only the final energy as output
    # furthermore, the number of iterations limits the optimization, thus we iterate over iterations, and 
    # save the energy ~ iterations
    bt <- TRUE
    if (bt){
      backtracking_iterations <- 10^(c(5, 10, 15)) # many zeros here :-) 
      
      RES_BACKTRACKING_TMP <- list()
      for (BT_ITER in 1:length(backtracking_iterations)){
        res_backtracking <- spectre:::run_optimization_backtracking(alpha_list = alpha_list,
                                                                    total_gamma = total_gamma[TOTAL_GAMMA],
                                                                    target = target_matrix_sim,
                                                                    max_runs = backtracking_iterations[BT_ITER],
                                                                    energy_threshold = 0.1)
        
        RES_BACKTRACKING_TMP[[BT_ITER]] <- cbind(n_sites[N_SITES],
                                      total_gamma[TOTAL_GAMMA],
                                      mean_alpha[TOTAL_GAMMA],
                                      res_backtracking$energy,
                                      backtracking_iterations[BT_ITER])
        
      }
      RES_BACKTRACKING[[TOTAL_GAMMA]] <- RES_BACKTRACKING_TMP
    }
  } # end of TOTAL_GAMMA-loop 
  # RES_MIN_CONF <-  RES_MIN_CONF 
  RES_MIN_CONF_ALL[[N_SITES]] <- RES_MIN_CONF 
  RES_RND_SWITCH_ALL[[N_SITES]] <- RES_RND_SWITCH
  RES_BACKTRACKING_ALL[[N_SITES]] <- RES_BACKTRACKING
} # end of N_SITES-loop

save(RES_MIN_CONF_ALL, file = "./data/RES_MIN_CONF_ALL.rda")
save(RES_RND_SWITCH_ALL, file = "./data/RES_RND_SWITCH_ALL.rda")
save(RES_BACKTRACKING_ALL, file = "./data/RES_BACKTRACKING_ALL.rda")
