### Plot results of algorithm tests...
library(ggplot2)

# random switching
# ================
load("./data/RES_RND_SWITCH_ALL.rda")
rnd_unlist <- unlist(RES_RND_SWITCH_ALL, recursive = FALSE, use.names = FALSE)
rnd_rbind <- data.frame(do.call(rbind, rnd_unlist))
names(rnd_rbind) <- c("sites", "gamma", "species", "energy", "iteration")

rnd_rbind$species <- rnd_rbind$species / rnd_rbind$gamma
rnd_rbind$gamma  <- as.factor(rnd_rbind$gamma)
rnd_rbind$species  <- as.factor(rnd_rbind$species)

ggplot(data = rnd_rbind, aes(x = iteration, y = energy, color = species))+
   geom_point(size = .2)+
  geom_smooth(method="loess", se=FALSE,span=0.3)+
  facet_grid(scales="fixed", vars(gamma),vars(sites),labeller = label_both)+
  labs(title = "random switching")+
  scale_x_continuous(name ="iterations",
                     breaks = seq(0, 5000, 2500))+
  scale_y_continuous(name ="energy",
                     limits = c(0,1),
                     breaks = seq(0, 1, 0.5))+
  theme_bw()

# "min_conf"
# ==========
load("./data/RES_MIN_CONF_ALL.rda")
min_unlist <- unlist(RES_MIN_CONF_ALL, recursive = FALSE, use.names = FALSE)
min_rbind <- data.frame(do.call(rbind, min_unlist))
names(min_rbind) <- c("sites", "gamma", "species", "energy", "iteration")

min_rbind$species <- min_rbind$species / min_rbind$gamma
min_rbind$gamma  <- as.factor(min_rbind$gamma)
min_rbind$species  <- as.factor(min_rbind$species)

ggplot(data = min_rbind, aes(x = iteration, y = energy, color = species))+
  geom_point(size = .2)+
  geom_smooth(method="loess", se=FALSE,span=0.3)+
  facet_grid(scales="fixed", vars(gamma),vars(sites),labeller = label_both)+
  labs(title = "min-conf")+
  scale_x_continuous(name ="iterations",
                     breaks = seq(0, 5000, 1000))+
  scale_y_continuous(name ="energy",
                     limits = c(0,1),
                     breaks = seq(0, 1, 0.5))+
  theme_bw()

# "backtracking"
# ==========
load("./data/RES_BACKTRACKING_ALL.rda")
back_unlist <- unlist(RES_BACKTRACKING_ALL, recursive = FALSE, use.names = FALSE)
back_unlist_2 <- unlist(back_unlist, recursive = FALSE, use.names = FALSE)
back_rbind <- data.frame(do.call(rbind, back_unlist_2))
names(back_rbind) <- c("sites", "gamma", "species", "energy", "iteration")

back_rbind$species <- back_rbind$species / back_rbind$gamma
back_rbind$gamma  <- as.factor(back_rbind$gamma)
back_rbind$species  <- as.factor(back_rbind$species)

ggplot(data = back_rbind, aes(x = log(iteration, base = 10), y = energy, color = species))+
  geom_point(size = 1.2)+
  geom_smooth(method="loess", se=FALSE,span=0.3)+
  facet_grid(scales="fixed", vars(gamma),vars(sites),labeller = label_both)+
  labs(title = "backtracking")+
  scale_x_continuous(name ="log(iterations, base=10)",
                     # breaks = seq(0, 5000, 100)
                     )+
  scale_y_continuous(name ="energy",
                     limits = c(0,max(back_rbind$energy)),
                     breaks = seq(0, 2, 0.5))+
  theme_bw()


