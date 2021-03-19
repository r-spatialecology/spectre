### Bray-Curtis vs. Commonness Optimization

library(ggplot2)
library(tidyverse)


### Bray-Curtis [branch=bc]

df<- list.files(path = paste0(getwd(), "/test_bc_MSP/res_bc"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_obs))+ 
  scale_y_continuous(name = expression(RCE(bc)) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)

ggplot(data = df) + 
  geom_boxplot(aes(y = bray_error))+ 
  scale_y_continuous(name = expression(BC_error(bc)) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)

ggsave(paste0(getwd(), "/test_bc_MSP/bc.png"), width = 8, height = 6, units = "in")


### Commonness [branch = dev_bc] 

df<- list.files(path = paste0(getwd(), "/test_bc_MSP/res_dev_bc"), pattern = "*.rds", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

df$replicate <- as.factor(df$replicate)

ggplot(data = df) + 
  geom_boxplot(aes(y = RCE_sol_obs))+ 
   scale_y_continuous(name = expression(RCE) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)

ggplot(data = df) + 
  geom_boxplot(aes(y = bray_error))+ 
  scale_y_continuous(name = expression(BC_error(dev_bc)) ) + #, breaks = seq(3.9, 4.0, 0.02 ), limits = c(3.91, 3.98)) + 
  scale_x_discrete(name ="" )+
  theme_bw(base_size = 18)

ggsave(paste0(getwd(), "/test_bc_MSP/dev_bc.png"), width = 8, height = 6, units = "in")


