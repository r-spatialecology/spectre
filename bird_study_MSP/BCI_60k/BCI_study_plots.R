### Plot results of min_conf_ algorithm tests...
setwd("~/spectre/bird_study_MSP/BCI_60k")
library(ggplot2)
library(tidyverse)

df<- list.files(path = getwd(), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$species <- as.factor(df$species)
df$replicate <- as.factor(df$replicate)

ggplot(data = df, aes(x = iteration, y = energy_before ))+ # , col = replicate))+
  geom_point(size = .3, alpha = 0.4)+
  facet_grid(scales="fixed", vars(species),vars(sites),labeller = label_both)+
  labs(title = "BCI, min-conf_0")+
  scale_x_continuous(name ="iterations", breaks = seq(0, 200000, 20000)) +
  scale_y_continuous(name = bquote("Discrepancy" ~D), breaks = seq(0, 1, 0.1), limits = c(0, .2) ) +
  theme_bw()

ggsave("BCI60kDisc.png", path = "../../vignettes", width = 8, height = 6, units = "in")

ggplot(data = df, aes(x = species, y = correct_pairs))+
  geom_point()+
  facet_grid(scales="fixed", vars(sites),labeller = label_both)+
  labs(title = "BCI , min-conf_0")+
  # scale_x_continuous(name ="iterations") +
  scale_y_continuous(name ="Correct pairs [%]", breaks = seq(0, 1, 0.5), limits = c(0, 1) ) +
  theme_bw()

ggsave("BCI60kPcor.png", path = "../../vignettes", width = 8, height = 6, units = "in")

# df$correct_pairs

