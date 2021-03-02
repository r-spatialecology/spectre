# Plot cumulated error in commonness 

library(ggplot2)
library(tidyverse)

df<- list.files(path = paste0(getwd(), "/results/BEN_2016_36sites"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

# boxplot RCE of solution (between solution and observed data) 

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_obs))+ 
  geom_hline(yintercept = df$RCE_target_observed[1], 
             linetype = "dotted", 
             size = 1.6, 
             color = "blue" ) + 
  scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)


ggsave("./figures/BEN_2016_sol_obs_boxpl.png", width = 8, height = 6, units = "in")


# boxplot RCE of solution (between solution and target):  RCE_sol_pred = (solution ~ predicted commonness)

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_pred))+ 
#  geom_hline(yintercept = df$RCE_target_observed[1], 
 #            linetype = "dotted", 
  #           size = 1.6, 
   #          color = "blue" ) + 
  scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)


ggsave("./figures/BEN_2016_sol_target_boxpl.png", width = 8, height = 6, units = "in")

# Evaluate prediction error(s)
# Errors introduced by the div models: (target - observed commonness)
print(paste0("The error between best estimate objective and observed commonness matrix was: ", df$RCE_target_observed[1]))

# Mean prediction error caused by the spectre algorithm
mean(df$RCE_sol_pred)

# Mean overall prediction error 
mean(df$RCE_sol_obs)


df_temp <- dplyr::group_by(df, replicate) %>%
  dplyr::summarize(obs_sol = mean(MAE_sol_obs),
                   tar_sol = mean(MAE_sol_pred))

# is the solution always closer to the observed data than the predicted target?
ggplot(data = df_temp) + 
  geom_boxplot(aes(y = obs_sol))+ 
  geom_hline(yintercept = df$MAE_target_observed[1], 
             linetype = "dotted", 
             size = 1.6, 
             color = "blue" ) + 
  scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)


ggsave("./manuscript/pics/allSolutionsBoxpl.png", width = 8, height = 6, units = "in")

### Error quantification

# target - solution (mean & sd)
round(mean(df_temp$tar_sol), 4)
round(sd(df_temp$tar_sol), 3)

# observations - solution (mean & sd)
round(mean(df_temp$obs_sol), 3)
round(sd(df_temp$obs_sol), 3)

## additional targets:
# pre_alpha + observed beta
df$MAE_pred_alpha_obs_beta[1]
# observed alpha and predicted beta
df$MAE_obs_alpha_pred_beta[1]

df$Commonness <- as.numeric(as.character(df$Commonness))

ggplot(data = df, aes(x = Commonness, y = Count/6.30, col = replicate))+ # , col = replicate))+
  facet_wrap(~replicate)  + 
  geom_point()+
  scale_y_continuous(name = "SiteXsite pairs [%]", breaks = seq(0, 100, 20), limits = c(0, 100) ) + 
  scale_x_continuous(name ="Commonness error", breaks = seq(0, 30, 1), limits = c(0, max(df$Commonness)))+
  theme_bw()

### keep only best solution

error_tmp <- unique(df$MAE_sol_pred)
best_replicate <- which.min(error_tmp)
print(paste0("Replicate ", best_replicate, " has the smallest (solution - gdm_target) error: " , round( error_tmp[best_replicate], 4)))

# create subset with best solution only
df <- df[df$replicate == best_replicate, ]

# overall error (observed - solution)
round(df$MAE_sol_obs[1], 3)

# target - solution
round(df$MAE_target_observed[1], 3)

print(paste0("The mean real commonness error was: ", round(df$MAE_sol_obs[1], 4)))
# print(paste0("The mean real sorensen error was: ", round(df$mean_real_sorensen_error[1], 4), 
#             ", the (solution - gdm_target) mean sorensen error was: ", round(df$mean_gdm_sorensen_error[1], 4) ))


# df$Commonness <- as.factor(df$Commonness)
ggplot(data = df, aes(x = Commonness, y = Count/6.30))+ # , col = replicate))+
  geom_point(size = 3)+
  scale_y_continuous(name = "SiteXsite pairs [%]", breaks = seq(0, 100, 20), limits = c(0, 100) ) + 
  scale_x_continuous(name ="Commonness error", breaks = seq(0, 30, 1), limits = c(0, max(df$Commonness)))+
  theme_bw(base_size = 18)


df$Count_cumsum <- cumsum(df$Count)
ggplot(data = df, aes(x = Commonness, y = Count_cumsum/6.30))+ # , col = replicate))+
  geom_point(size = 3)+
  scale_y_continuous(name = "SiteXsite pairs [%]", breaks = seq(0, 100, 20), limits = c(0, 100) ) + 
  scale_x_continuous(name ="Cumulative commonness error", breaks = seq(0, 100, 1), limits = c(0, max(df$Commonness)))+
  geom_vline(xintercept = df$MAE_sol_obs, linetype = "dotted", size = 1.6, color = "blue" ) +
  theme_bw(base_size = 18)
ggsave("../../manuscript/pics/MAE_sol_obs.png", width = 8, height = 6, units = "in")

df$MAE_target_observed
df$MAE_sol_obs

df$MAE_sol_pred

