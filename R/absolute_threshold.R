########################## absolute_threshold() ######################
#
#
#' Function to apply a given absolute threshold to the \code{.ncol} input graph
#' 
#' \code{absolute_threshold()} performs absolute thresholding on an input graph and 
#' writes the thresholded graph to an output file. This function is tailored specifically
#' for a clear absolute thresholding use case as compared to \code{threshold()}, which is the
#' more general purpose thresholding function. 
#' 
#' @param infile The \code{.ncol} graph file to be thresholded. 
#' @param outfile The output file path for the thresholded graph to be written to. 
#' @param threshold The desired absolute threshold to apply to the input graph.
#' @param overwrite A boolean parameter meant to prevent overwriting existing thresholded
#' graph files. The default (\code{FALSE}) will display a menu to the user if the passed
#' \code{outfile} already exists. Choosing \code{TRUE} will bypass this menu and overwrite
#' the existing file without interruption from a workflow.
#' @param sort_output Sorts the edges of the thresholded graph in descending order before
#' writing them to \code{outfile}. 
#' @examples
#' infile <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
#' thresholding::absolute_threshold(infile, 
#'                                  outfile = "./HCCS-ABSTHRESH.ncol",
#'                                  threshold = 0.90,
#'                                  sort_output = TRUE
#'                                  )
#' @returns Nothing. The thresholded graph is written to the file specified by outfile.
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
    
    t <- paste0("You are about to overwrite the graph file \'",
                outfile,
                "\'. Continue with absolute thresholding?")
    response <- menu(c("Yes", "No"), title = t)
                  
    # User does not want to overwrite output file, so end function   
    if(response == 2){
      print(paste0("--absolute_threshold() will not overwrite", outfile))
      print(paste0("--Leaving absolute_threshold() early..."))
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
