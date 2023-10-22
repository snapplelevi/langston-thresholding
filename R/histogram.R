#' edge_hist
#' 
#' Histogram function for displaying edge weight frequencies
#' 
#' @param infile tab separated .ncol formatted graph file (list of weighted edges - wel)
#' @param bin_width width of each histogram bin (default is 0.01)
#' @param sep determines how the wel file is separated (default is tab separated file)
#' @export
edge_hist <- function(infile,
                      bin_width=0.01,
                      sep = "\t"
                      )
{
  if(!file.exists(infile)){
    stop(paste("the file \"", infile, "\" does not exist"))
  }
  
  print(paste("Reading in the weighted edge list in ", 
              infile, "..."))
  
  # Read in data - handle errors?
  raw_data = read.table(infile, sep = sep, header=FALSE)
  print("Success - weighted edge list read!")
  
  # Error check if file is in proper format
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
  
  # Grab quantiles for histogram coloring
  quantiles <- quantile(raw_data$WEIGHT, 
                        prob=c(0.25,0.5,0.75), 
                        names=FALSE)
  
  ggplot2::ggplot(raw_data, ggplot2::aes(WEIGHT)) + 
    ggplot2::geom_histogram(binwidth=bin_width) +
    ggplot2::ggtitle(paste("Edge Weight Histogram for ", infile)) +
    ggplot2::xlab("Edge Weight") +
    ggplot2::ylab("Edge Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}