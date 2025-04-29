## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(thresholding)
input_graph_path <- system.file("extdata", "HumanCellCycleSubset.ncol", package="thresholding")

## -----------------------------------------------------------------------------
edge_hist_plot <- edge_hist(input_graph_path, title="Demo of edge_hist")
show(edge_hist_plot)

## -----------------------------------------------------------------------------
methods <- c(8, 3, 2)       # select desired analysis methods
starting_thresh <- 0.7   # choose lower bound thresholding value that thresholding loop 
num_samples <- 13        # running sig. and power - need the number of samples from data set

# Using overwrite = TRUE for vignette building purposes, read ?analysis for more 
# on this optional parameter
analysis_out_prefix <-analysis(input_graph_path, 
                               methods = methods, 
                               lower = starting_thresh,
                               num_samples = num_samples,
                               overwrite = TRUE
                               )

## -----------------------------------------------------------------------------
get_results_out <- get_results(analysis_out_prefix)
print(get_results_out)

## -----------------------------------------------------------------------------
chosen_threshold <- get_results_out[['scale_free']]
# or
chosen_threshold <- get_results_out$scale_free

print(chosen_threshold)

## -----------------------------------------------------------------------------
print(get_results_out[['TypeI-0.01']])
print(get_results_out[['Power-0.8']])

## ---- fig.width=7.5, fig.height=6---------------------------------------------
results_plot <- plot_t_vs_ev(analysis_out_prefix)

results_plot <- results_plot + ggplot2::ggtitle("Edge v. Vertices of Demo Graph")
show(results_plot)

## -----------------------------------------------------------------------------
# Save the output file path for convenience
output_graph_path = "./scale_free_threshold_0.82.ncol"
threshold(input_graph_path,
          outfile = output_graph_path,
          method = "absolute",
          thresh = chosen_threshold,
          overwrite = TRUE   # included to ensure graph file gets written in vignette
          )


## -----------------------------------------------------------------------------
# Display the results of absolute thresholding at 0.82
first_50_thresholded_lines <- utils::read.csv(output_graph_path, nrows=50)
print(first_50_thresholded_lines)

## -----------------------------------------------------------------------------
abs_output_graph_path = "./abs_thresh_demo_0.85"
absolute_threshold(input_graph_path,
                   abs_output_graph_path,
                   threshold = 0.85,
                   overwrite = TRUE,
                   )

# Display the results of absolute thresholding at 0.82
first_50_thresholded_lines <- utils::read.csv(abs_output_graph_path, nrows=50)
print(first_50_thresholded_lines)

