#' @title spectre
#'
#' @description
#' The goal of `spectre` is to provide an open source tool capable of predicting regional community composition at fine spatial resolutions using only sparse biological and environmental data.
#'
#' @name spectre
#' @docType package
#' @useDynLib spectre, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

globalVariables(c("error",
                  "i",
                  "value",
                  "x", 
                  "y"))
