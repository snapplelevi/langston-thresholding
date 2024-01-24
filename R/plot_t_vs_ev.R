#' Plots the edge (E) and vertex (V)count of a graph at varying thresholds
#' Marks non-null methods and their optimal thresholds against the 
#' V/E line plots
#'
#' @param plot_df Dataframe with threshold values - returned from 
#' output of get_iterative_results, which contains the detailed
#' analysis of the graph at each increment of the thresholding
#' process
#' @param D List of methods and their optimal thresholds from 
#' calling get_results. The user can either pass the resulting variable
#' from get_results or the list itself (i.e. D=variable$D instead of D=variable).
#' @export 
plot_t_vs_ev <- function(plot_df, D){
  
}