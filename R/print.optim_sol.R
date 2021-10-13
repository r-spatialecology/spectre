#' print.optim_sol
#'
#' @description Print method for optim_sol object
#'
#' @param x optim_sol object with randomized patterns.
#' @param ... Arguments passed to \code{cat}.
#'
#' @details
#' Printing method for optimized estimate of community composition.
#'
#' @seealso
#' \code{\link{run_optimization_min_conf}}
#'
#' @return void
#'
#' @aliases print.optim_sol
#' @rdname print.optim_sol

#' @export
print.optim_sol <- function(x, ...) {
  
  # get lowest error and total iterations
  lowest_error <- min(x$error$error)
  
  total_iterations <- max(x$error$i)
  
  # print result
  cat(paste0("Total iterations   : ", total_iterations, "\n",
             "Min absolute error : ", lowest_error, "\n",
             "Access 'optimized_grid' and 'error' for results."), ...)
}
