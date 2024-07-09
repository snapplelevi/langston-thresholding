# Function index (in order of appearance)
# 1. find_last()                    -> helper function , NOT exported
# 2. get_iterative_t_vals()         -> NOT exported
# 3. get_significance_t_vals()      -> NOT exported
# 4. get_local_global_alpha_value() -> NOT exported
# 5. get_results()                  -> exported

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
  
  
  print(files)
  # print(paste0("files inside: ", files))
  # Create array of data frames read in from files array
  all_dfs <- c()
  for(file in files){
    
    # Attempt to read in the file, but catch any errors and let user know
    # that the file was 
    tryCatch(
      
      expr = {
        
        df <- utils::read.csv(file, sep="")
        
      },
      error = function(err){
        writeLines("--get_results()  ---internal---> get_iterative_t_values():",  con = stderr())
        stop(paste0("the file: ", file, "\n
                     was not able to be read in by utils::read.csv.\n
                     Make sure the file path is correct, and there is valid
                     data in this file."))
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
    print('afterfffffff')
  }
  
  # Return empty df if there wasn't one big data frame made
  if(length(all_dfs) == 0){
    return(data.frame())
  }
  
  # Grouping by threshold and removing any duplicates from files with 
  # identical methods performed to preserve the max values
  df <- all_dfs %>% 
          dplyr::group_by(threshold) %>% 
          dplyr::filter(threshold == max(threshold)) %>% 
          dplyr::distinct()  %>% 
          dplyr::arrange(threshold)
  
  # print(colnames(df))
  
  writeLines("")
  writeLines("")
  # Reset index
  row.names(df) <- NULL

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
  diffs_b <- (abs(diffs) > drop_cutoff) & (diffs < 0)
  if(sum(diffs_b) > 0){
    # Find index of maximum difference and add the value of the
    # first index value (works same as df.index[0] in pandas)
    offset <- which.max(diffs_b) + strtoi(rownames(df)[1])
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
            dplyr::mutate(maximal.clique.count = 
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
    sub.app <- base::append(subdf$maximal.clique.count, NaN)
    
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
  # ------------OLD -------------
  # print("Nsv:")
  # print(Nsv)
  # print(sort(Nsv))
  # print("df$threshold")
  # print(df$threshold)
  # 
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
    if(sum(diffs > 0) > 3){
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
    
    if(!all(is.nan(df$random.clustering.coefficient))){
      # elo clustering coefficient
      # first local maximum
      #D['elo_clustering'] = df["threshold"][argrelextrema(C0_diffs, np.greater_equal)[0][0] + df.index[0]].min()
      
      tmp_diff <- df$clustering.coefficient - df$random.clustering.coefficient
      C0_diffs <- stats::runmed(tmp_diff, k = 3)
      elo_clustering_coefficient_diffs <- C0_diffs 
      df <- base::cbind(df, elo_clustering_coefficient_diffs)
      
      
      # Logic from python script (argrelextrema from scipy):
      #D['elo_clustering'] = df["threshold"][argrelextrema(C0_diffs, np.greater_equal)[0][0] + df.index[0]].min()
      #
      # Adapted to findpeaks, which has a different output logic. Please leave a Github issue
      # if this does not yield proper results. 
      peaks <- pracma::findpeaks(C0_diffs)
      D$D['elo_clustering'] <- min(df$threshold[ peaks[1][2] ] )
    }
    
  } # end outer if
  
  # Loop through D and only save the non NaN values (the methods the user requested
  # will have valid thresholds attributed to them
  # TODO: add this logic into Power/Sig and local_global (make function at top 
  # of this file perhaps)
  i = length(D$D)
  while(i > 0){
    if(is.nan(D$D[[i]])){
      D$D[[i]] <- NULL
    }
    i = i - 1
  }
  
  writeLines("############# get_iterative_t_values - DONE #############")
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
    
    line1 <- readLines(file, n=1)
    vals <- stringr::str_extract_all(line1, "\\d*(\\.)?\\d+")[[1]]
    
    alpha <- as.double(vals[1])
    sample_size <- as.integer(vals[[2]])
    r <- as.double(vals[[3]])
    
    D$D[paste0("TypeI-", alpha)] <- as.numeric(r)
    
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
  
  writeLines("############# get_significance_t_values - DONE #############\n")
  return(power_df)
}

################################################################################
#                    get_local_global_alpha_value() -- internal helper function
# Helper function for thresholding::get_results()
# Gets the alpha value from the local/global results
#   D_local_global is a part of Carissa's function implementation, but get_results
#   does not use it. For now, any interaction with 
get_local_global_alpha_value <- function(files, D_local_global=NULL){
  
  # Make empty list
  all_dfs <- c()
  
  # Append non-empty data frames to list
  for(file in files){
    
    df <- utils::read.csv(file, sep="", row.names=NULL) 
    
    # Attempt to read in the file, but catch any errors and let user know
    # that the file was 
    tryCatch(
      
      expr = {
        
        # row.names = 1 SAME index_col = 0 in python?
        df <- utils::read.csv(file, sep="", row.names=NULL)
        
      },
      error = function(err){
        writeLines("--get_results()  ---internal--->  get_local_global_alpha_value():", con = stderr())
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
      all_dfs <- rbind(all_dfs, df)
    }
  }
  
  if(length(all_dfs) == 0){
    return(data.frame())
  }
  
  #print(colnames(all_dfs))
  # Combine dfs and group by the alpha value
  # Equivalent of following line
  # power_df = pd.concat(all_power_df).groupby("r").max()#skipna=True)    
  df <- all_dfs %>%
        dplyr::group_by(alpha) %>%
        dplyr::filter(alpha == max(alpha)) %>%
        dplyr::distinct()
  
  # Reset index
  row.names(df) <- NULL
  
  # max_ac <- df[df["X2nd.eigenvalue"] < 1 & 
  #               df["almost.disconnected.component.count"] > 1]
  # max_ac <- subset(df, 
  #                  (X2nd.eigenvalue < 1) & (almost.disconnected.component.count > 1)
  #                  )
  #print(df$X2nd.eigenvalue)
  
  
  max_ac <- df %>% filter(X2nd.eigenvalue < 1 & almost.disconnected.component.count > 1)
  #print(max_ac)
  
  if(nrow(max_ac) > 0){
    min_alpha <- min(max_ac$alpha)
    row_alpha_exist <- subset(max_ac, alpha == min_alpha)

    alm_dis_max <- max(max_ac$almost.disconnected.component.count)
    
    r_a_m_tmp <- subset(max_ac, almost.disconnected.component.count == alm_dis_max)
    row_alpha_max <- subset(r_a_m_tmp, alpha == min(alpha))


    # Optional parameter D_local_global - add later but not necessary in core
    # functionality of get_results()
    # D_local_global["alpha_max"] = row_alpha_max.squeeze().to_dict()
    # D_local_global["alpha_exist"] = row_alpha_exist.squeeze().to_dict()
    
    # print(paste0("Row alpha exist: ", row_alpha_exist))
    # print(paste0("Row alpha max: ", row_alpha_max))
  }
  else{
    return(data.frame())
  }
  
  writeLines("############# get_local_global_alpha_value - DONE #############\n")
  return(df)
}

##############################################################################
#' Returns the resulting analysis method thresholding values
#' after running thresholding::analysis()
#' 
#' \code{get_results()} will take an output file prefix and attempt to accumulate
#' the results connected to that prefix. These results come in the form of analysis methods
#' and their recommended thresholding value (as determined by \code{get_results()}).
#' 
#' For example, take the output file prefix of "\code{EXAMPLE_NAME}". If this is passed
#' to \code{analysis()} once with \code{methods=c(4,5)}, and another time with
#' \code{methods=c(6,7)}, then there will be two output files:
#' \enumerate{
#'    \item EXAMPLE_NAME-45.iterative.txt
#'    \item EXAMPLE_NAME-67.iterative.txt
#' }
#' 
#' When given the prefix of "\code{EXAMPLE_NAME}", \code{get_results()} will look at these 
#' files and combine the results from each method into one data structure.
#' 
#' @param outfile_prefix File prefix for resulting output file from 
#' running the \code{analysis()} function (file would be <prefix>.iterative.txt)
#' @param plot_iterative Optionally plot the vertices and edges vs. threshold value as a graph.
#' uses ggplot2 to automatically call this package's plot_t_vs_ev() function without
#' the user needing to manually extract the required parameters.
#' 
#' @returns The returned list that can contain:
#' \itemize{
#'    \item A nested list of keyed on analysis method names. The values of these
#'    keys will be the recommended threshold from the corresponding method. This will 
#'    be called \code{"D"}.
#'    \item An alpha value representing the recommended significance value. This will 
#'    be called \code{"alpha"}
#' }
#' The output of \code{get_results()} and its separate wrappers will depend on the methods 
#' passed to \code{analysis()}. Values will either be valid or +/-\code{Inf}. If a method and/or value
#' for alpha does not show up in the returned list, then \code{get_results()} could not 
#' find this method or alpha value in the files with the matching prefix.
#' not show up 
#' @export
get_results <- function(outfile_prefix, plot_iterative = FALSE){
 
  # Make so that all file with outfile_prefix are fetched
  it_fnames <- Sys.glob(file.path(getwd(), 
                       paste0(outfile_prefix, "*.iterative.txt")))
  sig_fnames = Sys.glob(file.path(getwd(),
                         paste0(outfile_prefix, "*.statistical_errors.txt")))           
  loc_glo_fnames = Sys.glob(file.path(getwd(), 
                        paste0(outfile_prefix, "*.local_global.txt")))
  
  # Create list for each result type for convenient looping
  iterative <- list("iterative_result", it_fnames)
  significant <- list("significance_result", sig_fnames)
  local_global <- list("local_global_result", loc_glo_fnames)
  
  # Create list of lists
  all_results <- list(iterative, significant, local_global)
  
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
      # for(file in files){
      #   print(paste0("File #", i, ": ", file))
      #   i <- i + 1
      # }
      print("-------- Starting get_iterative_t_values --------")
      # Supress min() and max() warnings returning Inf
      df <- supWarn(get_iterative_t_values(files, D))
    } 
    else if(method == "significance_result"){
      print("-------- Starting get_significance_t_values --------")
      power_df <- supWarn(get_significance_t_values(files, D, min_power=0.8))
    } 
    else if(method == "local_global_result"){
      print("-------- Starting get_local_global_alpha value --------")
      df_and_alpha <- supWarn(get_local_global_alpha_value(files))
    }
  }
  
  writeLines("############# get_result - DONE #############\n")
  
  
  # Plot vertex and edge counts by threshold value if user specifies
  # instead of making seperate call to plot_t_vs_ev()
  if(plot_iterative == TRUE){
    writeLines("############# plot_t_vs_ev() called #############\n")
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
  return(list(D = D, alpha = alpha))
}


