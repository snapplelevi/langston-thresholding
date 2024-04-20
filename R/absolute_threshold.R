# Function to apply a given absolute threshold to the .ncol input graph
# Input:  .ncol file (same used in analysis function)
# Output: .ncol file with absolute threshold applied
#' @export
absolute_threshold <- function(infile,
                               outfile, 
                               threshold,
                               overwrite = FALSE,
                               sort_output = FALSE){
  
  
  # Ensure the provided file to be thresholded exists
  if(file.exists(infile) == FALSE){
    message(paste0("Could not find the file \'", infile, "\'."))
    stop("input file path was not found.")
  }
  
  
  # Present user with the option to overwrite the output file
  # if the overwrite parameter is left as FALSE
  # Set the parameter to TRUE for unconditional overwriting
  if(overwrite == FALSE && file.exists(outfile) == TRUE){
    
    t <- paste0("You are about to overwrite the output file \'",
                outfile,
                "\'. Continue with absolute thresholding?")
    response <- menu(c("Yes", "No"), title = t)
                  
    # User does not want to overwrite output file, so end function   
    if(response == 2){
      return(invisible(NULL))
    }
  }
  
  # Parameter checking based on threshold
  # Like eventually throw an error if value is too large / small or something
  
  # Read in input file / ensure that infile exists/accessible
  # Allow any form of white space as separator
  in_df <- read.table(infile, sep="", col.names=c("V1", "V2", "TH"))
  
  if(length(colnames(in_df)) != 3){
    message("The input file does not have three columns specified by the ncol format.")
    message("Proper file format should contain three columns like so:")
    message("      v1   v2  edge_weight")
    stop(paste0("Incorrect number of rows from input file (", infile, ")."))
  }
  # Choose rows within the bounds of the absolute threshold
  down_sel <- subset(in_df, (TH >= threshold) | (TH <= -1*threshold))
  
  # Sort the output from most negative weight to most positive weight
  # in case user wants to find most strongly connected vertices by looking
  # at the top or bottom of the file
  # Will sort in DESCENDING order by default
  if(sort_output){
    down_sel$TH <- base::sort(down_sel$TH, decreasing = TRUE)
  }
  
  # Write to the path specifed by outfile
  write.table(down_sel, 
              outfile, 
              row.names=FALSE,
              col.names=FALSE)
}

absolute_threshold("example/HumanCellCycleSubset.ncol", 
                   "example/TESTOKOKOK",
                   threshold=0.9,
                   overwrite=TRUE,
                   sort_output=TRUE
                   )