#' @title Thresholding Analysis
#' Thresholding functions - documentation handled through roxygen2 style comments
#'
#' @description
#' thresh.analysis - Primary Thresholding Analysis Function
#'
#' TODO: ADD PARAM DOCUMENTATION / REMOVE FROM FUNCTION PARAMS
#' @param outfile_prefix Prefix to outfile (Ex. "test." turns to "test.iterative.txt")
#'
#' @export
thresh.analysis <- function(outfile_prefix = "",
                            G,              # graph read from non-dimacs .wel graph file (igraph internally)
                            l = 0.5,        # Lower bound for thresholding
                            u = 0.99,       # Upper bound for thresholding
                            increment = 0.01,  # Threshold increment
                            windowSize = 5,    # Sliding window size for spectral method
                            minimumpartitionsize = 10,
                            minimum_cliquesize = 5,
                            min_alpha = 0,
                            max_alpha = 4,
                            alpha_increment = 0.1,
                            num_samples = 0,
                            significance_alpha = 0.01,
                            bonferroni_corrected = 0,
                            methods = ""      # Set of flagged methods
                            ){
    # Rest of function below
    print("Hello, graph thresholding! testing build")
    print(outfile_prefix)

    # Argument Checking?

    #### Calling the C++ function from the driven main analysis function
    # .Call("R_thresh_analysis", ...)
    #### Remember to cast parameters as a data type
    #### Pass on the graph file to the C code... let the error checking happen there


    # Print out Environment stuff?
}

#! FUTURE FUNCTIONS BELOW
