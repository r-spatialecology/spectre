#' @title plot_error
#' 
#' @description Plot the absolute error
#' 
#' @param x Results object from run_optimization_min_conf
#' 
#' @details 
#' Plot energy over time
#' 
#' @return ggplot
#' @references xxx

#' @export
plot_error <- function(x) {
  error_df <- x$error

  plot <- ggplot2::ggplot(data = error_df) + 
    ggplot2::geom_line(ggplot2::aes(x = i, y = error)) + 
    ggplot2::labs(x = "Iterations", y = "Error") +
    ggplot2::theme_bw()
  
  return(plot)
}