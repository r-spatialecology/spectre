### Plot results of min_conf_ algorithm tests...

library(ggplot2)
library(tidyverse)

df<- list.files(path = "./results", pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$species <- as.factor(df$species)

tabu_activated <- FALSE # only if tabu was "activated" in the test run 
if (tabu_activated){
  df$tabu <- df$tabu / df$sites 
  table(df$tabu)
  head(df)
  # filter by tabu value
  df <- filter(df, tabu == 0.2)
  names(df)
  table(df$tabu)
}

ggplot(data = df, aes(x = iteration, y = energy_scaled, color = species))+
  geom_point(size = .3, alpha = 0.4)+
  #geom_smooth(method="loess", se=FALSE,span=0.3)+
  facet_grid(scales="fixed", vars(gamma),vars(sites),labeller = label_both)+
  labs(title = "min-conf_1, tabu = 0%")+
  scale_x_continuous(name ="iterations") +
  scale_y_continuous(name ="energy", breaks = seq(0, 1, 0.2), limits = c(0, 1) ) +
  theme_bw()

ggplot(data = df, aes(x = iteration, y = energy_before, color = species))+
  geom_point(size = .3, alpha = 0.4)+
  #geom_smooth(method="loess", se=FALSE,span=0.3)+
  facet_grid(scales="fixed", vars(gamma),vars(sites),labeller = label_both)+
  labs(title = "min-conf, tabu = 20%")+
  scale_x_continuous(name ="iterations") +
  scale_y_continuous(name ="energy") +
  theme_bw()
