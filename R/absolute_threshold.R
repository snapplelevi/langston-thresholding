# Function to apply a given absolute threshold to the .ncol input graph
# Input:  .ncol file (same used in analysis function)
# Output: .ncol file with absolute threshold applied
#' @export
absolute_threshold <- function(infile, threshold, outfile=""){
  
  # Make sure output file doesn't already exist
  # Ask user if they would like to overwrite this file
  # if it does?
  # t <- paste0("You are about to overwrite output file \'", outfile,
  #                "\'. Continue?")
  # response <- menu(c("Yes", "No"), title=t)
  # if(response == nuh uh ) { leave the program}
  # Read in input file / ensure that infile exists/accessible
  
  
  # Parameter checking based on threshold
  
  
  # Results are stored in output file
  # Strip infile name and use that as default
  # Otherwise, just use the passed outfile name
}