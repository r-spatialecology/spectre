#' creates a dataframe of simulated environmental variable values for each landuse.
#' This is done by simple sampling from the data.
#' 
#' @param lutdf A table of site locations and landuses.
#' @param env A table of environmental variable values linked to each land use.
#' @return A dataframe of simulated environmental variable values for each landuse.

simulate_environmental_data <- function(lutdf, env) {
  gdmdat <- NULL
  # Loop trough env variables:
  for (i in (names(env)[-1])) {
    
    gdmdat.i <- NULL
    
    for (j in base::unique(env$PlotID))
    {
      ## Filter env:
      ## Extract the data for the current variable:
      env.ij <- env %>% dplyr::filter(PlotID == j) 
      env.ij <- env.ij[[i]]
      ## Split up lutdf
      lutdf.ij <- lutdf %>% dplyr::filter(lut == j) %>% dplyr::select(x, y)
      ## Create samples:
      lutdf.ij <- lutdf.ij %>% dplyr::mutate(z = sample(env.ij, nrow(lutdf.ij), replace=TRUE))
      ## Update names:
      names(lutdf.ij) <- c("x", "y", i)
      ##Attach to output:
      gdmdat.i <- rbind(gdmdat.i, lutdf.ij)
      
    }
    
    if (is.null(gdmdat)) {
      gdmdat <- gdmdat.i
    } else {
      ## Attach to output by joining:
      gdmdat <- gdmdat %>% dplyr::right_join(gdmdat.i, by=c("x", "y"))  
    }
  }
  
  ## Add cellid to envdat:
  gdmdat$cellid <- as.numeric(rownames(gdmdat))
  
  return(gdmdat)
}