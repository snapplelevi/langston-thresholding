#' edge_hist
#' 
#' Histogram function for displaying edge weight frequencies. \code{edge_hist()} requires that an 
  #' \link{https://lgl.sourceforge.net/} file be passed. In other words, the .ncol file can be thought
#' of as a \strong{w}eighted \strong{e}dge \strong{l}ist (\strong{wel}). The first column contains rows of vertices, the
#' second column contains rows of adjacent vertices, and the third column represents the weight for an edge
#' between the vertex in the first column one and the vertex in the second column. 
#' 
#' @param infile File path to the .ncol formatted graph file (list of weighted edges - \strong{wel})
#' @param bin_width (Optional) The width of each histogram bin. The default is bins of \strong{width} = 0.01
#' @param sep (Optional) Specifically decide how the columns in the \strong{wel} file are separated.
#' The default separator is any white space.
#' 
#' @examples
#' your_file_name_here <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
#' save_plot <- thresholding::edge_hist(your_file_name_here, bin_width = 0.01, sep = '\t')
#' show(save_plot)    ### the above call to edge_hist could also be wrapped in show()
#' 
#' @returns A ggplot object that can be displayed with show(), or by calling the edge_hist() directly in the R terminal.
#' @export
edge_hist <- function(infile,
                      bin_width=0.01,
                      sep = ""
                      )
{
  
  # Make sure the file exists 
  if(!file.exists(infile)){
    stop(paste("\redge_hist(): the file \"", infile, "\" does not exist.\nLeaving edge_hist()..."))
  }
  
  # Start message for user
  print(paste("Reading in the weighted edge list in ", 
              infile, "..."))
  
  # Read in .ncol file - FOR LATER: handle errors on 
  raw_data = read.table(infile, sep = sep, header=FALSE)
  print("Success - weighted edge list read!")
  
  # Error check if file is in proper format
  # This function assumes you pass in a three column .ncol file.
  if(ncol(raw_data) != 3){
    message("Error in edge_hist():")
    message(paste(infile, " is not in .ncol format. Expected threee columns (V1 V2 WEIGHT)"))
    message(paste0("Make sure the graph file is in .ncol format and that the sep parameter"))
    message(paste0("matches the separator used in the input file."))
    return(invisible(NULL))
  }
  
  
  # Rename the column names for standardizing the columns 
  colnames(raw_data)[1] <- "V1"
  colnames(raw_data)[2] <- "V2"
  colnames(raw_data)[3] <- "WEIGHT"
  
  # Validate that the third column is of type double or integer
  base::stopifnot(is.double(raw_data[["WEIGHT"]]) |
                  is.integer(raw_data[["WEIGHT"]]))
  
  #############   Future feature (beyond ver 1.0.0   ###############
  # Grab quantiles for histogram coloring
  # quantiles <- quantile(raw_data$WEIGHT, 
  #                       prob=c(0.25,0.5,0.75), 
  #                       names=FALSE)
  ##################################################################
  
  #############   MAKING THE PLOT   ###############
  # This returns the plot to the user to either
  # store in a variable or pass directly to a 
  # call to show()... ex.  show(edge_hist(infile))
  ggplot2::ggplot(raw_data, ggplot2::aes(WEIGHT)) + 
    ggplot2::geom_histogram(binwidth=bin_width,
                            fill="orange",
                            color="black") +
    ggplot2::ggtitle(paste("Edge Weight Histogram for ", infile)) +
    ggplot2::xlab("Edge Weight") +
    ggplot2::ylab("Edge Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}