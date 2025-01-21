# Function index (in order of appearance)
# 1. find_last()                    -> helper function , NOT exported
# 2. get_iterative_t_vals()         -> NOT exported
# 3. get_significance_t_vals()      -> NOT exported
# 4. get_results()                  -> exported

##############################################################################
#                         find_last()
# Helper function to find the index of the last occurrence of a charcter
# returns index of last character if found, returns -1 if not. 
#
# Used internally for finding multiple occurrences of file names the utility wrapper
# functions ( get_iter_t_vals(), get_sig_t_vals() )
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
#               get_iterative_t_values()  -- internal helper function
# Translated from Carissa Bleker's Thresholding helper functions
# for quicker analysis of thresholding analysis results
#
# The link to Carissa's Github repo: 
#   https://github.com/carissableker/thresholding
# The link to the combine_analysis_results functions: 
#   https://github.com/carissableker/thresholding/blob/master/example/combine_analysis_results.py

# Helper function for thresholding::get_results()
# Gets thresholding values from the iterative results
# Reads in tab the separated .iterative.txt files from 
# thresholding::analysis
# 
# @param files Vector of file paths or just one file path
# @param D Optional argument that accepts a list. This can be used
# to capture the thresholding results by a certain method. Not included in get_iter_t_vals()
# as of 2-27-24 for simplicity. Otherwise,get_results() fills out and returns this value automatically
# @param d_min_t INTERNAL CONTROL
get_iterative_t_values <- function(files,
                                   D=NULL,
                                   d_min_t=list(general=0)){
  
  

  # Create array of data frames read in from files array
  all_dfs <- c()
  for(file in files){
    
    # Attempt to read in the file, but catch any errors and let user know
    # that the file was not able to be opened
    tryCatch(
      
      expr = {
        
        df <- utils::read.csv(file, sep="")
        
      },
      error = function(err){
        writeLines("[get_results]:  get_results() ---internal---> get_iterative_t_values():",  con = stderr())
        stop(paste0("the file: ", file, "\n
                     was not able to be read in by utils::read.csv.\n
                     Make sure the file path is correct, and there is valid
                     data in this file. This can fail simply if the file is empty."))
      }
    )
    
    # if no rows, continue in loop
    writeLines("-------------------- Files and the number of rows -------------------- ")
    writeLines(paste0("File: ", file, "        nrows: ", nrow(df)))
    writeLines("")
    
    if(nrow(df) == 0)
      next
    
    # otherwise, append to array
    all_dfs <- rbind(all_dfs, df)
  }
  
  # Return empty df if there wasn't one big data frame made
  if(length(all_dfs) == 0){
    return(data.frame())
  }
  
  # Grouping by threshold and removing any duplicates from files with 
  # identical methods performed to preserve the max values
  # df <- all_dfs %>% 
  #         dplyr::group_by(threshold) %>% 
  #         dplyr::filter(threshold == max(threshold)) %>% 
  #         dplyr::distinct()  %>% 
  #         dplyr::arrange(threshold)
  
  # Group by "threshold" and select the maximum for each group
  df <- all_dfs %>%
        dplyr::group_by(threshold) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), max, na.rm = TRUE))
  
  # Sort by threshold
  df <- df %>% dplyr::arrange(threshold)
  
  # print(colnames(df))
  # print(df)
  
  # Reset index
  row.names(df) <- NULL
  
  # Formatting
  writeLines("")
  writeLines("")
  
  D$D['aplestin'] <- NaN
  D$D['cc_inflection'] <- NaN
  D$D['density_min'] <- NaN
  D$D['elo_clustering'] <- NaN
  D$D['gupta_clustering'] <- NaN
  D$D['mcr_2'] <- NaN
  D$D['mcr_3'] <- NaN
  D$D['mcr_max'] <- NaN
  D$D['rmt'] <- NaN
  D$D['scale_free'] <- NaN
  D$D['single_component'] <- NaN
  D$D['spectral_methods'] <- NaN
  D$D['whole_graph'] <- NaN
  
  if(nrow(df) < 3){
    return(df)
  }
  
  # Single component
  # First point before more than one cc appears
  D$D['single_component'] <- max(subset(df, 
                                  connected.component.count == 1,
                                  select = threshold
                                  ))
  
  # Whole graph
  D$D['whole_graph'] <- max(subset(df, 
                                vertex.count == max(vertex.count),
                                select = threshold
                                ))
  
  
  # Giant cc down inflection
  diffs <- diff(df$largest.cc.size)
  drop_cutoff <- sd(diffs)
  diffs_b <- ((abs(diffs) > drop_cutoff) & (diffs < 0))
  
  if(sum(diffs_b) > 0){
    # Find index of maximum difference and add the value of the
    # first index value (works same as df.index[0] in pandas)
    # Subtract 1 to account for R's 1-indexing...
    offset <- (which.max(diffs_b) - 1) + strtoi(rownames(df)[1])
    D$D['cc_inflection'] <- (df$threshold)[offset]
  }
  
  # density minimum
  D$D['density_min'] <- min(subset(df, 
                                 density == min(df$density),
                                 select = threshold))
  # scale free - lowest value where KS-pvalue > 0.1
  D$D['scale_free'] <- min(subset(df, 
                                scale.free.KS.p.value > 0.1,
                                select = threshold))
  
  # spectral methods: lowest t with largest number
  # of almost-disconnected-components
  # check if there are
  if("spectral" %in% names(d_min_t)){
    subdf <- subset(df, threshold > d_min_t['spectral'])
    # Reset index
    row.names(subdf) <- NULL
  }
  else{
    ########################################################
    # May need to change this if df is passed by reference #
    ########################################################
    subdf <- df
  }
              
  
  # Change values in subdf that are negative to be NaN
  # Otherwise, stay the same
  subdf <- subdf %>% 
           dplyr::mutate(almost.disconnected.component.count = 
                          dplyr::case_when(
                            almost.disconnected.component.count == -1 ~ NaN,
                            .default = almost.disconnected.component.count
                          ))
  
  # If there is at least one of the almost disconnected 
  # component counts that isn't NaN
  if(!all(is.nan(subdf$almost.disconnected.component.count))){
    # possiblilites: fiedler < 2 (?)
    # R prepends X since variables can't start with number
    max_ac <- subset(subdf, X2nd.eigenvalue < 2)
    if(nrow(max_ac) > 0){
      D$D['spectral_methods'] <- min(subset(max_ac, 
                                      almost.disconnected.component.count
                                      ==
                                      max(almost.disconnected.component.count),
                                      select = threshold))
    }
  }
  
  # RMT somewhere between 
  # poisson chi pvalue is > 0.05 (i.e. Poisson), 
  # goe chi2 pvalue > 0.05 (i.e. not GOE)
  if("rmt" %in% names(d_min_t)){
    subdf <- subset(df, threshold > d_min_t['rmt'])
    row.names(subdf) <- NULL
  }
  else{
    subdf <- df
  }
  
  
  poisson.pvalue <- min(subset(subdf, 
                               poisson.pvalue > 0.05, 
                               select=threshold))
  goe.pvalue <- max(subset(subdf, 
                           goe.pvalue < 0.05, 
                           select=threshold)) 
  D$D['rmt'] <- (poisson.pvalue + goe.pvalue) / 2 
  
  
  # Maximal clique ratio
  if("mcr" %in% names(d_min_t)){
    subdf <- subset(df, threshold > d_min_t['mcr'])
    row.names(subdf) <- NULL
  }
  else{
    subdf <- df
  }
  
  # Change values in subdf that are negative to be NaN
  # Otherwise, stay the same
  subdf <- subdf %>% 
            dplyr::mutate(maximal.clique.ratio = 
                            dplyr::case_when(
                              maximal.clique.count == -1 ~ NaN,
                              .default = maximal.clique.count
                            ))
  
  # Same logic as above but for maximal clique count
  if(!all(is.nan(subdf$maximal.clique.count))){
    ################################################
    #     Find out what this does in practice      #
    #     May need to tweak append logic           #
    ################################################
    sub.tmp <- subdf$maximal.clique.count
    sub.app <- base::append(subdf$maximal.clique.count[
                            2:length(subdf$maximal.clique.count)
                            ]
                            , NaN)
    
    base::transform(subdf,
                    maximal.clique.ratio = sub.tmp / sub.app)
    
    base::transform(subdf, 
                    maximal.clique.ratio = replace(
                      maximal.clique.ratio, Inf, NaN)
                    )

    D$D['mcr_2'] <- min(subset(subdf, maximal.clique.ratio >= 2,
                            select = threshold))
    
    D$D['mcr_3'] <- min(subset(subdf, maximal.clique.ratio >= 3,
                            select = threshold))
    
    D$D['mcr_max'] <- min(subset(subdf, 
                            maximal.clique.ratio 
                            == 
                            max(maximal.clique.ratio),
                            select = threshold))
  } # end if
  
  
  # aplestin
  # ------------OLD -------------  # 
  # Nsv <- df$edge.count / df$vertex.count
  
  # Error in pracma::gradient(Nsv, df$threshold) : 
  # Arguments 'h1' and 'h2' must be strictly increasing.
  #
  # dNsv_dt <- pracma::gradient(Nsv, df$threshold)
  # df <- base::cbind(df, dNsv_dt)
  
  
  # ----------- NEW ------------
  edge_unique <- base::unique(df$edge.count)
  vertex_unique <- base::unique(df$vertex.count)
  threshold_unique <- base::unique(df$threshold)
  
  Nsv <- edge_unique / vertex_unique
  Nsv_gradient <- pracma::gradient(Nsv, threshold_unique)
  dNsv_dt <- rep(Nsv_gradient, 
                  each = length(df$threshold) / length(threshold_unique)
                  )
    
  
  # Account for if the df has a different amount of ranges and repeated
  # thresholds (could happen through running different)
  rep_list <- (function(input_list){
      
      # Keep track of which indices in df$threshold are repeated so the same
      # index can be added N = number of repeats time to the returned vector
      dup_list <- duplicated(input_list)
      
      # First element of list will always be unique, so start at the second
      # for the while loop logic to be clean.
      rep_list <- c(1)
      ind <- 2
      inner_ind <- 2
      
      while(inner_ind <= length(dup_list)){
          
          streak <- 1
          
          while(inner_ind < length(dup_list) && dup_list[[inner_ind+1]] == TRUE){
              streak <- streak + 1
              inner_ind <- inner_ind + 1
          }
          
          rep_list <- append(rep_list, rep(ind, streak))
          ind <- ind + 1
          inner_ind <- inner_ind + 1
      }
      
      # Return list of repeated indices in passed vector (df$threshold)
      rep_list
  })(df$threshold)
  
  # Repeat the indices of the smaller nsv_gradient so that the same value 
  # gets assigned to each distinct threshold value
  matched_dNsv_dt <- Nsv_gradient[rep_list]
  
  # Insert new column into the data frame with the 'dNsv_dt' name
  # so an error doesn't pop up for the user
  df <- base::cbind(df, dNsv_dt = matched_dNsv_dt)
  D$D['aplestin'] <- min(df$threshold[dNsv_dt > 0])
  
  
  # same logic as above but for clustering coefficient
  if(!all(is.nan(df$clustering.coefficient))){
    # gupta clustering coefficient
    # largest t with sharp increase
    # estimate by first point where at least 
    # 3 differences are larger than 0.5 * stddev
    found_gupta <- FALSE
    diffs <- diff(df$clustering.coefficient)
    
    # first check: is there consistent increase in 
    # clustering coefficient?
    if(!any(is.na(diffs)) && sum(diffs > 0) > 3){
      cutoff <- sd(diffs[diffs > 0], na.rm=TRUE) * 0.3
      prev_prev_d <- diffs[1]
      prev_d <- diffs[2]
      i = 0
      for(d in 3:length(diffs)){
        values <- c(prev_prev_d, prev_d, diffs[[d]])
        if(all(values) > cutoff && i > 1){
          found_gupta <- TRUE
          break
        }
        prev_prev_d <- prev_d
        prev_d <- diffs[[d]]
        i = i + 1
      }
    } # end diffs check if
    
    if(found_gupta){
      D$D['gupta_clustering'] <- df$threshold[[i]]
    }

    
    # Make sure there are at least some non-nan or non-infinite values
    # in these clustering columns (this method was not requested if all
    # are NaN or Inf)
    if(!all(is.nan(df$random.clustering.coefficient)) 
       && !all(is.infinite(df$random.clustering.coefficient))
       && !all(is.nan(df$clustering.coefficient))
       && !all(is.infinite(df$clustering.coefficient))){
      # elo clustering coefficient
      # first local maximum
      
      #D['elo_clustering'] = df["threshold"][argrelextrema(C0_diffs, np.greater_equal)[0][0] + df.index[0]].min()
      
      tmp_diff <- df$clustering.coefficient - df$random.clustering.coefficient
      C0_diffs <- stats::runmed(tmp_diff, k = 3)

      elo_clustering_coefficient_diffs <- C0_diffs 
      df <- base::cbind(df, elo_clustering_coefficient_diffs)
      
      # Remove NaN values to prevent issues with pracma::findpeaks  
      C0_diffs <- C0_diffs[!is.nan(C0_diffs)]
      C0_diffs <- C0_diffs[!is.infinite(C0_diffs)]
      
      # Function to detect flat peaks
      find_flat_peaks <- function(v) {
        peaks <- pracma::findpeaks(v)
        
        # Initialize a list to store indices of flat peaks
        flat_peaks <- c()
        
        # Manually find the flat peaks
        for (i in 2:(length(v) - 1)) {
          if (v[i] == v[i - 1] || v[i] == v[i + 1]) {
            if (v[i] >= v[i - 1] && v[i] >= v[i + 1]) {
              flat_peaks <- c(flat_peaks, i)
            }
          }
        }
        
        return(list(peaks = peaks, flat_peaks = flat_peaks))
      }
      
      both_peaks <- find_flat_peaks(C0_diffs)
      is_peaks_null = is.null(both_peaks$peaks)
      is_flat_peaks_null = is.null(both_peaks$flat_peaks)
      
      
      offset <- NaN
      if(!is_peaks_null && !is_flat_peaks_null){
        offset <- min(both_peaks$flat_peaks[[1]], both_peaks$peaks[1][2])
      } else if (!is_peaks_null){
        offset <- both_peaks$peaks[1][2]
      } else if(!is_flat_peaks_null){
        offset <- both_peaks$flat_peaks[[1]]
      }
      
      D$D['elo_clustering'] = df$threshold[offset]
      
      # print(paste0("C0_diffs AFTER NaN removed:"))
      # print(C0_diffs)
      # 
      # # Adapted to findpeaks, which has a different output logic. Please leave a Github issue
      # # if this does not yield proper results
      # peaks <- pracma::findpeaks(C0_diffs)
      # print("peaks: ")
      # #print(peaks)
      # #print(paste0("min expression : ", min(df$threshold[ peaks[1][2] ] )))
      # D$D['elo_clustering'] <- min(df$threshold[ peaks[1][2] ] )
    }
    
  } # end outer if
  

  writeLines("[get_results]:  get_iterative_t_values() - DONE \n")
  #print(paste0("D after: ", D$D))
  
  return(df)
}


################################################################################
#                     get_significance_t_values()  -- internal helper function
# Helper function for thresholding::get_results()
# Gets the significance thresholding values from significance results
get_significance_t_values <- function(files, D, alpha=0.5, min_power=0.8){
 
  # Create list of all power dfs
  all_power_df <- c()  
  
  for(file in files){
    # Attempt to read in the file, but catch any errors and let user know
    # that the file was not able to be opened
    tryCatch(
      
      expr = {
        
        df <- utils::read.csv(file, sep="")
        
      },
      error = function(err){
        writeLines("[get_results]:  get_results() ---internal---> get_significance_t_values():",  con = stderr())
        stop(paste0("the file: ", file, "\n
                     was not able to be read in by utils::read.csv.\n
                     Make sure the file path is correct, and there is valid
                     data in this file. This can fail simply if the file is empty."))
      }
    )

    # Print info to user about file
    writeLines("-------------------- Files and the number of rows -------------------- ")
    writeLines(paste0("File: ", file, "        nrows: ", nrow(df)))
    writeLines("")

    line1 <- readLines(file, n=1)
    vals <- stringr::str_extract_all(line1, "\\d*(\\.)?\\d+")[[1]]
    
    alpha <- as.double(vals[1])
    sample_size <- as.integer(vals[[2]])
    r <- as.double(vals[[3]])
    
    D$D[paste0("TypeI-", alpha)] <- as.numeric(r)
    writeLines(paste0("\tType I: ", as.numeric(r), '\n'))
    # Attempt to read in the file, but catch any errors and let user know
    # that the file was 
    tryCatch(
      
      expr = {
        
        ####### figure out if index_col=0 is the default ######
        df <- utils::read.csv(file, sep="", skip=2)
        
      },
      error = function(err){
        writeLines("--get_results()  ---internal--->  get_significance_t_values():", con = stderr())
        stop(paste0("the file: ", file, "\n
                     was not able to be read in by utils::read.csv.\n
                     Make sure the file path is correct, and there is valid
                     data in this file."))
      }
    )
    
    if(nrow(df) == 0){
      next
    }
    else{
      all_power_df <- rbind(all_power_df, df)
    }
    
  } # end file for-loop
  
  if(length(all_power_df) == 0){
    D$D[paste0("TypeI-", alpha)] <- NaN
    D$D[paste0("Power-", min_power)] <- NaN
    return(data.frame())
  }
  
  
  # Equivalent of following line
  #power_df = pd.concat(all_power_df).groupby("r").max()#skipna=True)    
  power_df <- all_power_df %>%
    dplyr::group_by(r) %>%
    dplyr::filter(r == max(r)) %>%
    dplyr::distinct()
  
  row.names(power_df) <- NULL
  
  D$D[paste0("Power-", min_power)] <- min(power_df["r"][power_df["power"] >= min_power])
  
  
  writeLines("[get_results]: get_significance_t_values() - DONE \n")
  return(power_df)
}


##############################################################################
#' Returns the resulting analysis method thresholding values
#' after running \code{thresholding::analysis()}
#' 
#' @description
#' \strong{Note: modifying the column names of the output file(s) produced by} \code{analysis()} 
#' \strong{will lead to unintended behavior.}
#' 
#' \code{get_results()} will take an output file prefix and attempt to accumulate
#' the results connected to that prefix. This prefix come in the form of analysis methods
#' and their recommended thresholding value (as determined by \code{get_results()}).
#' 
#' @details
#' For example, take the output file prefix of "\code{EXAMPLE_NAME}". If this is passed
#' to \code{analysis()} once with \code{methods=c(4,5)}, and another time with
#' \code{methods=c(6,7)}, then there will be two output files:
#' \enumerate{
#'    \item EXAMPLE_NAME-45.iterative.txt
#'    \item EXAMPLE_NAME-67.iterative.txt
#' }
#' When given the prefix of "\code{EXAMPLE_NAME}", \code{get_results()} will look at these 
#' files and combine the results from each method into one data structure.
#' 
#' @param outfile_prefix string. Prefix for the resulting output file from 
#' running the \code{analysis()} function (file would be <prefix>.iterative.txt). The 
#' output files will be assumed to be in the current working directory from where
#' \code{get_results()} was called if no path to the prefix is specified. 
#' @param plot_iterative boolean. Optionally plot the vertices and edges vs. threshold value.
#' This uses ggplot2 to automatically call this package's \code{plot_t_vs_ev()} function without
#' the user needing to manually extract the required parameters.
#' @param return_dfs boolean. Returns the iterative and significance data frames if they exist. These 
#' are returned by \code{get_iter_t_vals()} and \code{get_sig_t_vals()} independently as well.
#' Defaults to \code{FALSE} since these data frames can get rather long. The outputted list will
#' be structured differently based on this input. 
#' 
#' @returns The returned list that can contain:
#' \itemize{
#'    \item A nested list of keyed on analysis method names. The values of these
#'    keys will be the recommended threshold from the corresponding method. "
#'    \item A nested list containing data frames and a list of methods list. Them methods list
#'    is the same as described above. The data frames are included if they are valid based on the
#'    \code{analysis()} output files found.
#' }
#' The output of \code{get_results()} and its separate wrappers will depend on the methods 
#' passed to \code{analysis()}. Values will be valid (non-\code{NaN} and non-\code{Inf}). 
#'
#' @examples
#' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
#' file.copy(data_file, "./")     # Copy the file to your working directory
#' outfile_prefix = "./get_results_test"
#' analysis(data_file, 
#'          outfile_prefix = outfile_prefix,
#'          methods = c(2,3,8),
#'          num_samples = 13,
#'          overwrite = TRUE
#'          )
#' thresholding::get_results(outfile_prefix, plot_iterative = TRUE)
#' @export
get_results <- function(outfile_prefix, 
                        plot_iterative = FALSE,
                        return_dfs = FALSE){
 
  # Make so that all file with outfile_prefix are fetched
  it_fnames <- Sys.glob(file.path(paste0(outfile_prefix, "*.iterative.txt")))
  sig_fnames = Sys.glob(file.path(paste0(outfile_prefix, "*.statistical_errors.txt")))           

  # Create list for each result type for convenient looping
  iterative <- list("iterative_result", it_fnames)
  significant <- list("significance_result", sig_fnames)
  
  # Create list of lists
  all_results <- list(iterative, significant)
  
  D <- new.env()  # Environment passed to helper functions (persists values)
  D$D <- list()   # Hash-map/Dictionary of method/threshold value pairs
  alpha <- NaN  # Alpha significance level
  
  # Used later to suppress warnings in analysis helper functions. These functions
  # use max() and min() from base, which return Inf/-Inf when the vector 
  # passed to them is empty (which may arise if a certain condition is not met)
  # though analysis. 
  #
  # We use a temporary wrapper workaround for the longe "suppressWarnings" by 
  # storing the function name in a shorter variable.
  supWarn <- suppressWarnings
  
  # Loop through all results
  for(res in all_results){
    method <- res[[1]]  # string value
    files <- res[[2]]   # list of strings
    
    if(method == "iterative_result"){
      i <- 1

      writeLines("[get_results]:  Starting get_iterative_t_values() \n")
      # Supress min() and max() warnings returning Inf
      df <- supWarn(get_iterative_t_values(files, D))
    } 
    else if(method == "significance_result"){
      writeLines("[get_results]:  Starting get_significance_t_values() \n")
      power_df <- supWarn(get_significance_t_values(files, D, min_power=0.8))
    } 
  }
  
  writeLines("[get_results]:  DONE  \n\n")
  
  # Loop through D and only save the non NaN values (the methods the user requested
  # will have valid thresholds attributed to them
  # TODO: add this logic into Power/Sig
  i = length(D$D)
  while(i > 0){
    if(is.nan(D$D[[i]]) || is.infinite(D$D[[i]])){
      D$D[[i]] <- NULL
    }
    i = i - 1
  }
  
  # Plot vertex and edge counts by threshold value if user specifies
  # instead of making seperate call to plot_t_vs_ev()
  if(plot_iterative == TRUE){
    writeLines("\n[get_results]: plot_t_vs_ev() called\n")
    show(thresholding::plot_t_vs_ev(outfile_prefix))
  }
  
  # Results accessible by <out_variable>$D ; <out_variable>$alpha
  # or
  #               <out_variable>["D"] ; <out_variable>["alpha"]
  #
  # ex. OUTPUT <- get_results("FILE-PREFIX")
  #     to get whole list:         OUTPUT$D
  #     to get a specific element: OUTPUT$D$cc_inflection
  
  D <- D$D
  
  if(return_dfs){
    r_list <- list()
    # Include iterative df if valid
    if(all(dim(df)) != 0){
      r_list$iter_df <- df
    }
    
    # Include significance df if valid
    if(all(dim(power_df)) != 0){
      r_list$sig_df <- power_df
    }
    # Include the methods list (regardless) 
    r_list$methods <- D
    
    return(r_list)
  }
  
  # Just return the methods list if return_dfs = F
  return(D)
}


