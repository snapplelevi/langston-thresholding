#' Prints the resulting analysis method thresholding values
#' after running thresholding::analysis()
#' 
#' Note: This function assumes that the .iterative anaylsis() output file
#' is in the current directory unless the path to the file is provided. 
#' 
#' @param outfile filename or file path for resulting output file from 
#' running the analysis function (file would be <prefix>.iterative.txt)
get_results <- function(outfile_prefix){
  # Make so that all file with outfile_prefix are fetched
  fname <- file.path(outfile)
  
  if(!file.exists(fname)){
    stop(paste("the file ", fname, " does not exist."))
  }
  
  
}

