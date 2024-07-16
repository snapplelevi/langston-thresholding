##############################################################################
#                         find_last()
# Helper function to find the index of the last occurrence of a charcter
# returns index of last character if found, returns -1 if not. 
#
# Used internally for finding multiple occurrences of file names the utility wrapper
# functions ( get_iter_t_vals(), get_sig_t_vals(), get_local_global_alpha() )
find_last <- function(str, str_to_find){
  l = nchar(str)
  ind = -1
  
  while(l > 0){
    
    if(substr(str, l, l) == str_to_find){
      ind = l
      break
    }
    
    l <- l - 1
  }
  
  return(ind)
}

##############################################################################
#'                          get_iter_t_vals()
#'                          
#' User wrapper function for get_iterative_t_values
#' 
#' Returns the thresholding data frame created by the internal
#' get_iterative_t_values, which does not return anything to the user. This
#' data frame contains step-by-step thresholding figures for analysis methods, 
#' edge counts, vertex counts, and other metrics computed as a result of the 
#' \code{analysis()} function.
#' 
#' 
#' 
#' @param outfile_prefix Prefix of output file, which can have several
#' output file paths if the same prefix is run with \code{analysis()} several times. 
#' 
#' Ex.) \code{get_iter_t_vals("iter-prefix")} will execute \code{get_iterative_t_values}
#' on all files names \code{iter-prefix-###.iterative.txt}, where \code{###} is the unique combination
#' of methods passed to \code{analysis()}. This value can be equal to nothing ("") if no
#' methods were passed to \code{analysis()}.
#' 
#' @param recursive Option to allow for subdirectory searching of the file prefix.
#' Set to FALSE for default in case of large subdirectories. This option is meant to help 
#' save time with R Studio and other editors' auto completion feature, which allows for tabbing
#' to complete the file name from the present working directory (PWD) and maintaining the current
#' PWD. 
#' 
#' If \code{FALSE}, \code{get_iter_t_vals()} will use the path passed to the function from outfile_prefix.
#' This may be either relative or absolute. Only the directory containing the desired file(s)
#' will be searched. (Note: if no file paths are explicitly provided, the PWD (".") will be used.)
#' 
#' If \code{TRUE}, \code{get_iter_t_vals()} will recursively searching through all files matching the Regex pattern
#' of \code{^output_prefix.*\\.iterative\\.txt$}, which essentially allows for any file name that starts 
#' with the specified outfile_prefix and ends with ".iterative.txt".
#' 
#' @returns A list containing both a dataframe and a list. The dataframe contains the raw interative
#' thresholding data at each threshold value for the requested methods. The nested list contains
#' each requested method (for the particular analysis output file(s) found from the output prefix)
#' and its corresponding threshold value.
#' @export
get_iter_t_vals <- function(outfile_prefix, recursive=FALSE){
  
  
  # Strip out the file path (if it exists) so the subsequent Regex
  # works properly
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
  patt <- paste0("^", outfile_prefix, ".*\\.iterative\\.txt$")
  
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
    stop("\rExiting get_iter_t_vals...")
  }
  
  # Run the iterative analysis on the found files
  D_iter <- new.env()
  D_iter$D <- list()
  iter_df <- suppressWarnings(get_iterative_t_values(it_fnames, D_iter))

  return(list(iter_df = iter_df, D_iter = D_iter))
}
