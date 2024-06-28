##############################################################################
#'                          get_sig_t_vals()
#' User wrapper function for \code{get_significance_t_values()}
#' 
#' Returns the thresholding data frame created by the internal
#' \code{get_significance_t_values()}, which does not return anything to the user.
#' 
#' The returned data frame includes graph and graph method values for each 
#' increment of the threshold. 
#' 
#' 
#' @param outfile_prefix Prefix of output file, which can have several
#' output file paths if the same prefix is run with analysis() several times. 
#' 
#' Ex.) get_iter_t_vals("iter-prefix") will execute get_iterative_t_values
#' on all files names iter-prefix-###.iterative.txt where ### is the process
#' ID of the internal process (like 4318, 3341, 414143, etc.).
#' 
#' @param recursive Option to allow for subdirectory searching of the file prefix.
#' Set to FALSE for default in case of large subdirectories. This option is meant to help 
#' save time with R Studio and other editors' auto completion feature, which allows for tabbing
#' to complete the file name from the present working directory (PWD) and maintaining the current
#' PWD. 
#' 
#' If \code{FALSE}, \code{get_sig_t_vals()} will use the path passed to the function from outfile_prefix.
#' This may be either relative or absolute. Only the directory containing the desired file(s)
#' will be searched. (Note: if no file paths are explicitly provided, the PWD (".") will be used.)
#' 
#' If \code{TRUE}, \code{get_sig_t_vals()} will recursively searching through all files matching the Regex pattern
#' of \code{^output_prefix.*\\.iterative\\.txt$}, which essentially allows for any file name that starts 
#' with the specified outfile_prefix and ends with ".iterative.txt".
#' @export
get_sig_t_vals <- function(outfile_prefix, recursive=FALSE){
  
  # Strip out the file path (if it exists) so the subsequent Regex
  # works properly
  # find_last defined in "get_results.R"
  path_end <- find_last(outfile_prefix, "/")
  path <- "."
  
  # Strip out directory path if there was a final '/' found in the outfile_prefix
  if(path_end > 0){
    
    if(recursive==FALSE){
      path <- base::substr(outfile_prefix, 1, path_end)
    }
    
    outfile_prefix <- base::substr(outfile_prefix, path_end+1, nchar(outfile_prefix))
    
  }
  
  # Create regex pattern to match files with exact .iterative.txt format
  # Allows for file names with prefix to be changed, but must keep the .iterative.txt
  # format to work in this case.
  patt <- paste0("^", outfile_prefix, ".*\\.statistical_errors\\.txt$")
  
  
  # Find all possible files with the given file prefix
  it_fnames <- list.files(path=path, 
                          recursive=recursive, 
                          pattern=patt,
                          full.names = TRUE)
  
  # No files found for the given prefix
  if(length(it_fnames) == 0){
    pwd_mess <- ifelse(recursive, "PWD and subdirectories", "PWD")
    message(paste0("No files found in ", pwd_mess, " for file prefix: \"", outfile_prefix, "\" ."))
    message("Please check the file prefix and try again.")
    message("Iterative files will have the format of: outfile_prefix-####.iterative.txt")
    message("   - Note: #### represents the integer value of the process ID of the")
    message("           running process when the iterative file was created from analysis().")
    message("")
    message("You may also try using the recursive=TRUE parameter in order to search recursively")
    message("through the current working directory's subdirectories for your file.")
    message("")
    message("If this doesn't work, please check the spelling of the file prefix.")
    message("")
    stop("\rExiting get_sig_t_vals")
  }
  
  # Run the iterative analysis on the found files
  sig_df <- suppressWarnings(get_significance_t_values(it_fnames))
  
  return(sig_df)
}