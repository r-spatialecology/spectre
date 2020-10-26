### Plot results of min_conf_ algorithm tests...

setwd("~/spectre/bird_study_MSP/BCI_known_sites")

library(ggplot2)
library(tidyverse)

df<- list.files(pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$species <- as.factor(df$species)
df$replicate <- as.factor(df$replicate)
df$known_sites <- as.factor(df$known_sites)
names(df)

ggplot(data = df, aes(x = iteration, y = energy_before, col = known_sites ))+ # , col = replicate))+
  geom_point(size = .4, alpha = 0.6)+
  facet_grid(scales="fixed", vars(species),vars(sites),labeller = label_both)+
  labs(title = "BCI, min-conf_0")+
  scale_x_continuous(name ="iterations", breaks = seq(0, 200000, 10000)) +
  scale_y_continuous(name = bquote("Discrepancy" ~D), breaks = seq(0, 1, 0.1), limits = c(0, .2) ) +
  theme_bw()

ggsave("BCI60kDiscKnownSites.png", path = "../../vignettes", width = 8, height = 6, units = "in")

ggplot(data = df, aes(x = known_sites, y = correct_pairs, col = species))+
  geom_point()+
  facet_grid(scales="fixed", vars(sites),labeller = label_both)+
  labs(title = "BCI , min-conf_0")+
  # scale_x_continuous(name ="iterations") +
  scale_y_continuous(name ="Correct siteXsite pairs [%]", breaks = seq(0, 100, 20), limits = c(50, 105) ) +
  theme_bw()

ggsave("BCI60kPcorKnownSites.png", path = "../../vignettes", width = 8, height = 6, units = "in")

# correct species ~ known sites
ggplot(data = df, aes(x = known_sites, y = correctly_predicted_species, col = species))+
  geom_point()+
  facet_grid(scales="fixed", vars(sites),labeller = label_both)+
  labs(title = "BCI , min-conf_0")+
  # scale_x_continuous(name ="iterations") +
  scale_y_continuous(name ="Correctly predicted species", breaks = seq(0, 1, .2), limits = c(0, 1) ) +
  theme_bw()

ggsave("BCI60kSpeciesKnownSites.png", path = "../../vignettes", width = 8, height = 6, units = "in")


