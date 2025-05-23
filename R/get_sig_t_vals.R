##############################################################################
#                          get_sig_t_vals()
#
#' User wrapper function for the internal \code{get_significance_t_values()}
#' 
#' @description
#' \strong{Note: modifying the column names of the output file(s) produced by} \code{analysis()} 
#' \strong{will lead to unintended behavior.}
#' 
#' Returns the thresholding data frame created by the internal
#' \code{get_significance_t_values()}.
#' 
#' @returns
#' The returned data frame includes Power analysis for values of 'r' and the corresponding
#' power. 
#' 
#' 
#' @param outfile_prefix string. Prefix of output file produce from \code{analysis()} 
#' (identical to this parameter in this package's \code{get_results()}), which can have several
#' output file paths if the same prefix is run with \code{analysis()} several times. 
#' 
#' Ex.) \code{get_sig_t_vals("sig-prefix")} will execute \code{get_significance_t_values()}
#' on all files names \code{sig-prefix-###.iterative.txt} where \code{###} is the unique combination
#' of methods passed to \code{analysis()}. This value can be equal to nothing ("") if no
#' methods were passed to \code{analysis()}.
#' 
#' @param recursive boolean. Option to allow for subdirectory searching of the file prefix.
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
#' 
#' @examples
#' library(ThresholdTuner)
#' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "ThresholdTuner") 
#' outfile_prefix <- tempfile('get_sig_t_vals')
#' 
#' analysis(data_file, 
#'          outfile_prefix = outfile_prefix,
#'          methods = c(2,3,8),
#'          num_samples = 13,
#'          overwrite = TRUE
#' )
#' get_sig_t_vals(outfile_prefix)
#' 
#' @export
get_sig_t_vals <- function(outfile_prefix, recursive=FALSE){
  # Find all possible files with the given file prefix (same as get_results)
  sig_fnames = Sys.glob(file.path(paste0(outfile_prefix, "*.statistical_errors.txt")))           
  
  # No files found for the given prefix
  if(length(sig_fnames) == 0){
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
  sig_df <- suppressWarnings(get_significance_t_values(sig_fnames, list()))
  
  return(sig_df)
}