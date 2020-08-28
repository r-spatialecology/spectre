### Plot results of min_conf_ algorithm tests...

setwd("~/spectre/BCI_tidy")

library(ggplot2)
library(tidyverse)

df<- list.files(path = paste0(getwd(), "/table_res"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

df$Commonness <- as.numeric(as.character(df$Commonness))

names(df)

df$delta_commonness <- df$mean_real_commonness_error / df$mean_commonness * 100

# MAE_commonness
ggplot(data = df, aes(x = iterations, y = mean_real_commonness_error, col = replicate))+ # , col = replicate))+
  facet_grid(scales="fixed", vars(n_species),vars(n_sites),labeller = label_both)+
  geom_point()+
  scale_y_continuous(name = "mean absolute commonness error") + # , breaks = seq(0, 100, 2), limits = c(0, 10) ) + 
  scale_x_continuous(name ="Max iterations", breaks = seq(0, 20000, 10000), limits = c(0, 20000))+
  theme_bw()

# delta_commonness error = MEA_commonness / mean absolute commonness 
ggplot(data = df, aes(x = iterations, y = delta_commonness, col = replicate))+ # , col = replicate))+
  facet_grid(scales="fixed", vars(n_species),vars(n_sites),labeller = label_both)+
  geom_point()+
  scale_y_continuous(name = "relative commonness error [%]") + # , breaks = seq(0, 100, 2), limits = c(0, 10) ) + 
  scale_x_continuous(name ="Max iterations", breaks = seq(0, 20000, 10000), limits = c(0, 20000))+
  theme_bw()

# distance D
ggplot(data = df, aes(x = iterations, y = final_distance_D, col = replicate))+ # , col = replicate))+
  facet_grid(scales="fixed", vars(n_species),vars(n_sites),labeller = label_both)+
  geom_point()+
  scale_y_continuous(name = "Distance D") + # , breaks = seq(0, 100, 2), limits = c(0, 10) ) + 
  scale_x_continuous(name ="Max iterations", breaks = seq(0, 20000, 10000), limits = c(0, 20000))+
  theme_bw()

names(df)
# distance D ~ mean_richness
ggplot(data = df, aes(x = mean_richness, y = final_distance_D, col = replicate))+ # , col = replicate))+
  facet_grid(scales="fixed", vars(iterations),labeller = label_both)+
  geom_point()+
  scale_y_continuous(name = "Distance D") + # , breaks = seq(0, 100, 2), limits = c(0, 10) ) + 
  scale_x_continuous(name ="Mean richness" ) + # , breaks = seq(0, 20000, 10000), limits = c(0, 20000))+
  theme_bw()
