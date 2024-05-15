#' edge_hist
#' 
#' Histogram function for displaying edge weight frequencies. `edge_hist()` requires that an 
#' [.ncol](https://lgl.sourceforge.net/) file be passed. In other words, the .ncol file can be thought
#' of as a \strong{w}eighted \strong{e}dge \strong{l}ist (\strong{wel}). The first column contains rows of vertices, the
#' second column contains rows of adjacent vertices, and the third column represents the weight for an edge
#' between the vertex in the first column one and the vertex in the second column. 
#' 
#' @param infile tab separated .ncol formatted graph file (list of weighted edges - \strong{wel})
#' @param bin_width width of each histogram bin (default is 0.01)
#' @param sep determines how the columns in the \strong{wel} file are separated (default is any white space)
#' 
#' @examples
#' your_file_name_here <- './example/HumanCellCycleSubset.ncol'
#' thresholding::edge_hist(your_file_name_here, bin_width = 0.01, sep = '\t')
#' 
#' @returns A ggplot object that can be displayed with show(), or by calling the function directly in the R terminal.
#' @export
edge_hist <- function(infile,
                      bin_width=0.01,
                      sep = ""
                      )
{
  
  # Make sure the file exists 
  if(!file.exists(infile)){
    stop(paste("the file \"", infile, "\" does not exist"))
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
    message(paste(infile, " is not in .ncol format. Expected threee columns (V1 V2 WEIGHT)"))
    stop(paste("the file \"", infile, "\" is not in .ncol format."))
  }
  
  
  # Rename the column names for standardizing the columns 
  colnames(raw_data)[1] <- "V1"
  colnames(raw_data)[2] <- "V2"
  colnames(raw_data)[3] <- "WEIGHT"
  
  # Validate that the third column is of type double or integer
  stopifnot(is.double(raw_data$WEIGHT) |
            is.integer(raw_data$WEIGHT))
  
  # Grab quantiles for histogram coloring - feature for later
  quantiles <- quantile(raw_data$WEIGHT, 
                        prob=c(0.25,0.5,0.75), 
                        names=FALSE)
  
  
  ############   MAKING THE PLOT   #############
  ggplot2::ggplot(raw_data, ggplot2::aes(WEIGHT)) + 
    ggplot2::geom_histogram(binwidth=bin_width,
                            fill="orange",
                            color="black") +
    ggplot2::ggtitle(paste("Edge Weight Histogram for ", infile)) +
    ggplot2::xlab("Edge Weight") +
    ggplot2::ylab("Edge Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}