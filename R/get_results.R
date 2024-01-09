# Translated from Carissa Bleker's Thresholding helper functions
# for quicker analysis of thresholding analysis results
#
# The link to Carissa's Github repo: 
#     https://github.com/carissableker/thresholding
# The link to the combine_analysis_results functions: 
#     https://github.com/carissableker/thresholding/blob/master/example/combine_analysis_results.py

# Helper function for thresholding::get_results()
# Gets thresholding values from the iterative results
# Reads in tab the separated .iterative.txt files from 
# thresholding::analysis
get_iterative_t_values <- function(files,
                                   D,
                                   d_min_t=list(general=0)){
  print(paste0("D before: ", D$D))
  # print(paste0("files inside: ", files))
  # Create array of data frames read in from files array
  all_dfs <- c()
  for(file in files){
    df <- read.csv(file, sep="\t")
    
    # if no rows, continue in loop
    print(nrow(df))
    if(nrow(df) == 0)
      next
    # otherwise, append to array
    all_dfs <- rbind(all_dfs, df)
  }
  
  # Return empty df if there wasn't one big data frame made
  if(length(all_dfs) == 0){
    return(data.frame())
  }
  
  # Grouping by threshold and removing any duplicates
  # to preserve the max values
  df <- all_dfs %>% 
          dplyr::group_by(threshold) %>% 
          dplyr::filter(threshold == max(threshold)) %>% 
          dplyr::distinct()  %>% 
          dplyr::arrange(threshold)
  
  # Reset index
  row.names(df) <- NULL
  print(colnames(df))
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
  D$D['Power-0.8'] <- NaN        
  D$D['TypeI-0.01'] <- NaN
  
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
    #######################################################
    # May need to change this if df is passed by reference#
    #######################################################
    subdf <- df
  }
                
  print(subdf$almost.disconnected.component.count)
  
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
                           select=threshold)) / 2
  D$D['rmt'] <- poisson.pvalue + goe.pvalue
  
  
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
  Nsv <- df$edge.count / df$vertex.count
  dNsv_dt <- pracma::gradient(Nsv, df$threshold)
  df <- base::cbind(df, dNsv_dt)
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
  
  print("##### get_iterative_t_values - DONE #####")
  #print(paste0("D after: ", D$D))
  
  return(df)
}



# Helper function for thresholding::get_results()
# Gets the significance thresholding values from significance results
get_significance_t_values <- function(files, D, alpha=0.5, min_power=0.8){
 
  # Create list of all power dfs
  all_power_df <- c()  
  
  for(file in files){
    
    lines <- readLines(file)
    line1 <- lines[[1]]
    vals <- stringr::extract_all(line1, "\\d*(\\.)?\\d+")
    
    alpha <- vals[[1]]
    sample_size <- vals[[2]]
    r <- vals[[3]]
    D$D[paste0("TypeI-", alpha)] <- as.numeric(r)
    
    ####### figure out if index_col=0 is the default ######
    df <- read.csv(file, sep="\t", skip=2, row.names=1)
    
    if(nrow(df) == 0){
      next
    }
    else{
      all_power_df <- append(all_power_df, df)
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
  
  D$D[paste0("Power-", min_power)] <- min(power_df[power_df["power"] >= min_power]["r"])
  
  print("##### get_significance_t_values - DONE #####")
  return(power_df)
}



# Helper function for thresholding::get_results()
# Gets the alpha value from the local/global results
#   D_local_global is a part of Carissa's function implementation, but get_results
#   does not use it. For now, any interaction with 
get_local_global_alpha_value <- function(files, D_local_global=NULL){
  
  # Make empty list
  all_dfs <- c()
  
  # Append non-empty data frames to list
  for(file in files){
    # row.names = 1 SAME index_col = 0 in python?
    df <- read.csv(file, sep="\t", row.names=1) 
    
    if(nrow(df) == 0){
      next
    }
    else{
      all_dfs <- append(all_dfs, df)
    }
  }
  
  if(length(all_dfs) == 0){
    return(data.frame())
  }
  
  # Combine dfs and group by the alpha value
  # Equivalent of following line
  # power_df = pd.concat(all_power_df).groupby("r").max()#skipna=True)    
  df <- all_dfs %>%
        dplyr::group_by(alpha) %>%
        dplyr::filter(alpha = max(alpha)) %>%
        dplr::distinct()
  
  # Reset index
  row.names(df) <- NULL
  
  max_ac <- df[df["X2nd-eigenvalue"] < 1 & 
                df["almost-disconnected-component-count"] > 1]
  max_ac <- subset(df, 
                   (X2nd-eigenvalue < 1) & (almost-disconnected-component-count > 1)
                   )
  
  if(nrow(max_ac) > 0){
    min_alpha <- min(max_ac$alpha)
    row_alpha_exist <- max_ac[max_ac$alpha == min_alpha]
    
    
    alm_dis_max <- max(max_ac$almost.disconnected.component.count)
    r_a_m_tmp <- max_ac[max_ac$almost.disconnectd.component.count == alm_dis_max]
    row_alpha_max <- r_a_m_tmp[r_a_m_tmp$alpha == min(r_a_m_tmp$alpha)]
    

    # Optional parameter D_local_global - add later but not necessary in core
    # functionality of get_results()
    # D_local_global["alpha_max"] = row_alpha_max.squeeze().to_dict()
    # D_local_global["alpha_exist"] = row_alpha_exist.squeeze().to_dict()
    
    print(paste0("Row alpha exist: ", row_alpha_exist))
    print(paste0("Row alpha max: ", row_alpha_max))
  }
  else{
    return(data.frame())
  }
  
  print("##### get_local_global_alpha_value - DONE")
  return(df)
}



#' Prints the resulting analysis method thresholding values
#' after running thresholding::analysis()
#' 
#' Note: This function assumes that the .iterative anaylsis() output file
#' is in the current directory unless the path to the file is provided. 
#' 
#' @param outfile_prefix filename or file path for resulting output file from 
#' running the analysis function (file would be <prefix>.iterative.txt)
#' @export
get_results <- function(outfile_prefix){
 
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
  
  # Loop through all results
  for(res in all_results){
    method <- res[[1]]  # string value
    files <- res[[2]]   # list of strings
    
    if(method == "iterative_result"){
      i <- 1
      for(file in files){
        print(paste0("File #", i, file))
        i = i + 1
      }
      print("Running get_iterative_t_values")
      df <- get_iterative_t_values(files, D)
    } 
    else if(method == "significance_result"){
      print("Running get_significance_t_values")
      power_df <- get_significance_t_values(files, D, min_power=0.8)
    } 
    else if(method == "local_global_result"){
      print("Running get_local_global_alpha value")
      df_and_alpha <- get_local_global_alpha_value(files)
    }
  }
  
  print("get_result - DONE")
  
  # accessible by <out_variable>$D ; <out_variable>$alpha
  # or
  #               <out_variable>["D"] ; <out_variable>["alpha"]
  #
  # ex. OUTPUT <- get_results("FILE-PREFIX")
  #     whole list:       OUTPUT$D
  #     specific element: OUTPUT$D$cc_inflection 
  D <- D$D
  return(list(D = D, alpha = alpha))
}
