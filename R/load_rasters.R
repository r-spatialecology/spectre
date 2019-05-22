#' Function to help find and load required raster files
#' 
#' @param path The directory path to the required data.
#' @param pattern The type of file the rasters are stored as.
#' @return Raster stack of imported files.
#' 
#' @import raster

load_raster <- function(path, pattern) {
  files <- list.files(path, pattern=pattern)
  maps <- raster::stack(paste(path,files, sep = ""))
  return(maps)
}