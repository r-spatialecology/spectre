# Plot cumulated error in commonness 

library(ggplot2)
library(tidyverse)

### spectre dev 

df<- list.files(path = paste0(getwd(), "/test_bc_MSP/res"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_pred))+ 
  #  geom_hline(yintercept = df$RCE_target_observed[1], 
  #            linetype = "dotted", 
  #           size = 1.6, 
  #          color = "blue" ) + 
  scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)


ggsave(paste0(getwd(), "/test_bc_MSP/dev_.png"), width = 8, height = 6, units = "in")

### Bray curtis 

df<- list.files(path = paste0(getwd(), "/test_bc_MSP/res_bc"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_pred))+ 
  #  geom_hline(yintercept = df$RCE_target_observed[1], 
  #            linetype = "dotted", 
  #           size = 1.6, 
  #          color = "blue" ) + 
  scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)


ggsave(paste0(getwd(), "/test_bc_MSP/bc_.png"), width = 8, height = 6, units = "in")


