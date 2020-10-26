### if empirical data is used in spectre, accuracy of solutions probably depends
### both on the spectre-algorithm, as well as on the structure of the data set.
### Thus this script visualizes/summarizes the data

### get data from betapart package
library("betapart") 
data(bbsData)
data <- bbs2000 # either 1980 or 2000 (years)

data[1:10, 1:10]
# data is formatted as: sites/states (rows); species (columns) 

dim(data) # 49 sites, 569 species 

# richness per site
(alpha_list <- rowSums(data))
# gamma diversity
(total_gamma <- sum(colSums(data) > 0))
plot( 1:length(alpha_list), alpha_list, ylim = c(0, total_gamma + 1),
      xlab = "Site", ylab = "Richness", main = "Blue = gamma-diversity, red = mean(richness)")
abline(h = total_gamma, lty = "dotted", lwd = 2, col = "blue")
abline(h = mean(alpha_list), col = "red", lty = 2)

### mean_richness / gamma ratio
mean(alpha_list) # 2000: 177
sd(alpha_list) # 2000: 43 

mean(alpha_list) / total_gamma # 0.31 




# abline(h = median(alpha_list), col = "green")


hist(alpha_list, breaks = 20, main = "Richness histogramm", xlab = "Richness")

### commonness matrix
commonness_explore <- as.matrix(spectre:::calculate_solution_commonness_rcpp(data))

plot_commonness_explore <- function(commonness_explore) {
  
  # commonness_explore <- calculate_solution_commonness_rcpp(species_grid[[1]])
  
  n_row <- nrow(commonness_explore)
  n_col <- ncol(commonness_explore)
  
  #commonness_explore[upper.tri(commonness_explore, diag = TRUE)] <- NA
  
  commonness_difference_grid <- commonness_explore 
  
  commonness_difference <- expand.grid(x = seq_len(n_row), 
                                       y = seq_len(n_col))
  
  commonness_difference$value <- commonness_difference_grid[seq_len(n_row * n_col)]
  
  plot <- ggplot2::ggplot(data = commonness_difference) + 
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value)) + 
    
    ggplot2::scale_fill_gradient2(name = "Commonness", 
                                  # midpoint = mean(max(commonness_difference$value)),
                                  low = "red", 
                                  mid = "white", 
                                  high = "blue") + 
    #scale_colour_gradient2(colours = rainbow(5)) + 
    ggplot2::coord_equal() + 
    ggplot2::theme_void()
  
  return(plot)
}

plot_commonness_explore(commonness_explore)


