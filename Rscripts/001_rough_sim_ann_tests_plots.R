in_dir <- "/results/001_SA001"
out_dir <- "/figures/001_SA001"

library("ggplot2")
library(tidyverse)

df<- list.files(path = paste0(getwd(), in_dir), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

names(df)

# RCE ~ T_0_factor
(plot1 <- ggplot(data = df, aes(x = factor(T_0_factor), y = RCE, col = factor(presence_prob)))+ # , col = replicate))+
    facet_grid(~presence_prob)+
    geom_boxplot()+
    ylab("RCE [%]")+
    theme_bw())
# ggsave(paste0(getwd(), out_dir, "a.png"))

(plot2 <- ggplot(data = df, aes(x = T_0_factor, y = min_error, col = factor(T_0_factor)))+ # , col = replicate))+
    facet_grid(~presence_prob)+
    geom_boxplot()+
    ylab("min_error")+
    theme(axis.text.x = element_text(color="#993333", 
                                     size=12, angle=90))) # +
    #theme_bw())

