# tests whether rescaled max and min values do not match the
# proportions of gamma entered
context("generate_data_raster_rescale")

testthat::test_that("Rescaled values have expected minimum and maximum values", {
  library(NLMR)
  library(raster)
  
  gamma <- 500
  min_proportion <- 0.1
  min_gamma <- gamma * min_proportion
  max_proportion <- 0.4
  max_gamma <- gamma * max_proportion
  
  test_landscape <- NLMR::nlm_random(ncol=50, 
                                     nrow=50)
  
  rescaled_landscape <- generate_data_raster_rescale(landscape, 
                                                     gamma, 
                                                     min_proportion, 
                                                     max_proportion)
  
  # Test min
  testthat::expect_equal(raster::minValue(rescaled_landscape),
                         min_gamma,
                         label = "Minimum rescaled gamma")
  
  # Test max
  testthat::expect_equal(raster::maxValue(rescaled_landscape),
                         max_gamma,
                         label = "Maximum rescaled gamma")
})
