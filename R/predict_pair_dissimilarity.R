#' Function creates a predicted pairwise dissimilarity for every pair of sites in the landscape
#' 
#' @param x x-coord.
#' @param y y-coord.
#' @param envdat A dataframe of environmental variable values for each landuse.
#' @param gdmmodel GDM model.
#' @return A pairwise matrix of predicted dissimilarity.

## This function creates a predicted pairwise dissimilarity for every pair of sites in the landscape

## envdat must be in the form created by simulate_environmental_variables.R
predict_pair_dissimilarity <- function(x, y, envdat, gdmmodel) 
{
  ## Extract data from envdat:
  #envdat.xy <- envdat %>% dplyr::filter(cellid %in% c(x, y)) %>% dplyr::select(-cellid)
  envdat.xy <- rbind(envdat[x,], envdat[y,]) %>% dplyr::select(-cellid)
  envdat.xy$plotID <- seq_len(nrow(envdat.xy))
  
  ### Select to front and remove lon and lat column:
  envdat.xy <- envdat.xy %>% dplyr::select(plotID, dplyr::everything()) %>% dplyr::select(-lon, -lat)
  envdat.m <- as.matrix(envdat.xy, dimnames=names(envdat))
  envdat.sp <- envdat.m[,1:4]
  gdm.tab <- gdm::formatsitepair(envdat.sp, bioFormat = 1, abundance=F, dist="bray", siteColumn = "plotID", XColumn = "x", YColumn = "y", predData = envdat.m)
  
  ## Do predictions:
  gdm.tab.pred <- gdm::predict.gdm(gdmmodel, gdm.tab, time=FALSE)
  
  return(gdm.tab.pred)  ## Write to matrix
} 