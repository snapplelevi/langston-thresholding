#' edge_hist
#' 
#' Histogram function for displaying edge weight frequencies
#' 
#' @param infile tab separated .ncol formatted graph file (list of weighted edges - wel)
#' @param bin_width width of each histogram bin 
#' @param sep determines how the wel file is separated (default is tab separated file)
#' @export
edge_hist <- function(infile,
                      bin_width,
                      sep = "\t"
                      )
{
  if(!file.exists(infile)){
    stop(paste("the file \"", infile, "\" does not exist"))
  }
  
  print(paste("Reading in the weighted edge list in ", 
              infile, "..."))
  
  # Read in data - handle errors?
  raw_data = read.table(infile, sep = sep, header=TRUE)
  print("Success - weighted edge list read!")
  
  # Print data
  print(raw_data)
}