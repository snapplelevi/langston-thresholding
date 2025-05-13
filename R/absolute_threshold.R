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
#' @param infile string. The path to the \code{.ncol} graph file to be thresholded. 
#' @param outfile string. The output file path for the thresholded graph to be written to. 
#' @param threshold numeric. The desired absolute threshold to apply to the input graph.
#' @param overwrite boolean. A parameter meant to prevent overwriting existing thresholded
#' graph files. The default (\code{FALSE}) will display a menu to the user if the passed
#' \code{outfile} already exists. Choosing \code{TRUE} will bypass this menu and overwrite
#' the existing file without interruption from a workflow.
#' @param sort_output boolean. Sorts the edges of the thresholded graph in descending order before
#' writing them to \code{outfile} if set to \code{TRUE}. Otherwise, descending order is not guaranteed. 
#' @examples
#' library(thresholding)
#' infile <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
#' outfile <- tempfile(fileext = ".ncol")
#' absolute_threshold(infile, 
#'                   outfile = outfile,
#'                   threshold = 0.90,
#'                   sort_output = TRUE,
#'                   overwrite = TRUE
#'                   )
#' print(paste0("Thresholded graph in '", outfile, ".'"))
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
    response <- utils::menu(c("Yes", "No"), title = t)
                  
    # User does not want to overwrite output file, so end function   
    if(response == 2){
      print(paste0("--absolute_threshold() will not overwrite", outfile))
      print(paste0("--Leaving absolute_threshold() early..."))
      return(invisible(NULL))
    } else {
      print(paste0("----Continuing with absolute_threshold()."))
      print(paste0("----Overwriting the graph file: ", outfile))
    }
  }
  
  # Parameter checking based on threshold
  # Like eventually throw an error if value is too large / small or something
  
  # Read in input file / ensure that infile exists/accessible
  # Allow any form of white space as separator
  vertex_col_1 <- "V1"
  vertex_col_2 <- "V2"
  threshold_col <- "threshold_column"
  in_df <- utils::read.table(infile, sep="", col.names=c(vertex_col_1,
                                                         vertex_col_2,
                                                         threshold_col))
  
  if(length(colnames(in_df)) != 3){
    message("The input file does not have three columns specified by the ncol format.")
    message("Proper file format should contain three columns like so:")
    message("      v1   v2  edge_weight")
    stop(paste0("Incorrect number of rows from input file (", infile, ")."))
  }
  # Choose rows within the bounds of the absolute threshold
  down_sel <- subset(in_df, 
                     (threshold_column >= threshold) | (threshold_column <= -1*threshold)
                     )
  
  # Sort the output from most negative weight to most positive weight
  # in case user wants to find most strongly connected vertices by looking
  # at the top or bottom of the file
  # Will sort in DESCENDING order by default
  if(sort_output){
    down_sel$threshold_column <- base::sort(down_sel$threshold_column, decreasing = TRUE)
  }
  
  # Write to the path specifed by outfile
  utils::write.table(down_sel, 
                     outfile, 
                     row.names=FALSE,
                     col.names=FALSE)
}
