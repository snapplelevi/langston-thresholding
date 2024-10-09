#' Plots the edge (E) and vertex (V) counts of a graph at varying thresholds
#' 
#' Plots a graph's edge and vertices by threshold value on two different scales. This function is intended to 
#' give more visual insight to the thresholds chosen by the iterative analysis methods of \code{analysis()},
#' which are methods 3 through 8. Non-NaN and non-Inf analysis methods are plotted at their optimal thresholds
#' against the V/E line plots. 
#' 
#' 
#' @details 
#' The vertical, orange lines represent a thresholding analysis method. The intersection point of both the V and E
#' thresholding curves are denoted by red and blue squares. The method name is annotated for each vertical line. The x value (s)
#' (threshold amount) that the method has yielded are where the vertical lines will appear. Use \code{get_results()} to see the 
#' concrete value of these thresholds.
#' 
#' @param iter_prefix string. The prefix of the .iterative.txt file(s) from the output \code{analysis()} function, which
#' is printed to the screen after \code{analysis()} ends.
#' @returns Nothing. Plots the ggplot2 object directly from the function.
#' @examples
#' library(thresholding)
#' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
#' file.copy(data_file, "./")     # Copy the file to your working directory
#' outfile_prefix = "./plot_t_vs_ev_test"
#' analysis(data_file, 
#'          outfile_prefix = outfile_prefix,
#'          methods = c(2,3,8),
#'          num_samples = 13,
#'          overwrite = TRUE
#'          )
#' plot_t_vs_ev(outfile_prefix)
#' @export 
plot_t_vs_ev <- function(iter_prefix){
  
  
  # # Validate that iter_prefix is valid
  # if(file.exists(iter_prefix) == FALSE){
  #   writeLines(paste0("-- plot_t_vs_ev(): The input file \'", iter_prefix, "\'"))
  #   writeLines(paste0("--                 was not able to be opened... please check spelling or the"))
  #   writeLines(paste0("--                 working directory if the file path is relative."))
  #   writeLines(paste0("--"))
  #   writeLines(paste0("-- Leaving plot_t_vs_ev()"))
  #   return(invisible(NULL))
  # }
  # 
  # 
  # # Strip the .iterative.txt suffix from the precise file name to use in get_iter_t_vals
  # iter_ext <- ".iterative.txt"
  # stripped_list <- unlist(strsplit(iter_prefix, iter_ext, fixed=TRUE))
  # 
  # # Make sure that there is only one .iterative.txt and that the extension
  # # is truly the extension (end of the file name)
  # if(length(stripped_list) != 2 &&
  #    stripped_list[[2]] != iter_ext){
  #   print("Print out error messages for .iterative.txt file ext. being wrong")
  #   return(invisible(NULL))
  # }
  # 
  #
  
  
  # Get the results from analysis to get the D$D list
  #stripped_prefix <- stripped_list[[2]]
  iter_results <- thresholding::get_iter_t_vals(iter_prefix)
  plot_df <- iter_results$iter_df
  
  
  # # Make sure get_results had something to plot
  # # Otherwise just tell the user there is nothing to plot
  # if(is.null(iter_results$D)){
  #   writeLines(paste0("-- plot_t_vs_ev(): internal use of get_results(\'", iter_prefix, "\') returned"))
  #   writeLines(paste0("--                 no valid analysis method thresholds. This means there"))
  #   writeLines(paste0("--                 is nothing to plot."))
  #   writeLines(paste0("--"))
  #   writeLines(paste0("-- Leaving plot_t_vs_ev()"))
  #   return(invisible(NULL))
  # }
  
  # Store D_iter$D list into D for easier reading / less typing
  D <- iter_results$D_iter$D
  
  # Set up structures for storing plotting data
  annotations <- list()
  labels <- list()
  names_vec <- names(D)
  
  ###########################################################
  # NEED TO REMOVE POWER AND SIGNFICANCE AND ANY OTHER 
  # ARTIFACT THINGS ---- THIS IS ONLY FOR ITERATIVE RESLUTS!!1
  # THE OTHER THINGS SCREW UP THE PLOTTING BIGGGG TIME
  
  # Used for keeping track of which element of the list you're on
  #     Used to get the string representation of the analysis method name
  #     for each corresponding value in D
  index <- 1
  
  # Loop through threshold values stored in D
  for(cur_val in D){
    
    # Extract the string 'key' for the current threshold val in D
    method_name <- names(D)[index]
    
    # Skip these two fields of D since they are non-iterative
    # They will screw up the plotting later if left in the labels
    if (grepl("^Power-", method_name)==TRUE|| grepl("TypeI-", method_name)==TRUE){
      next;
    }
    

    # If the method produced a non-zero value, add the method name (annotation)
    # and its corresponding value to the plotting sets.
    if(is.nan(cur_val) == FALSE && is.infinite(cur_val) == FALSE){
      
      is_in <- FALSE

      # Used to access annotation  name in annotations list
      annot_index <- 1

      for(annot_val in annotations){
        # String value of annotation name
        annot <- names(annotations)[annot_index]
        
        # START DEBUGGING IF THIS WORKS OR NOT
        if(method_name == 'spectral_methods'){
          
          # print(paste0(cur_val, " and ", annot_val))
          # print(paste0("method_name: ", method_name, " and ", annot))
          # print(dplyr::near(cur_val, annot_val, tol=1e-09))
        }        
        # END DEBUGGING BLOCK
        
        
        if(dplyr::near(cur_val, annot_val, tol=1e-09)){
          # print(paste0("adding ", method_name, " to the ", annot, " vector"))
          labels[[annot]] <- c(labels[[annot]], method_name)
          annotations[[annot]] <- mean(cur_val, annot_val)
          # Make sure the element does not get added to the labels
          # since it is close to the current value and already in the set
          is_in <- TRUE
        } 
        
        # Increment name index as well
        annot_index <- annot_index + 1
      }

      # Only add if the method hasn't been seen before
      if(is_in == FALSE){
        # print(paste0("adding: ", method_name))
        labels[[method_name]] <- c(method_name)   # start vector of labels, hashed on name
        annotations[[method_name]] <- cur_val    # value that corresponds to label name(s)
      }
      # print(paste0("labels after", method_name))
      # print(labels)
      # print(paste0("annotations after ", method_name))
      # print(annotations)
    }  # end outer is.nan if

    # Increment the index to get the string name of the method
    index <- index + 1

  } # end for loop
  
  # 
  # print("Annotations: ")
  # print(annotations)
  # 
  # print("labels: ")
  # print(labels)
  
  #print(head(plot_df))
  
  v_count <- plot_df$vertex.count[1]
  e_count <- plot_df$edge.count[1]
  t_begin <- plot_df$threshold[1]
  t_end <- utils::tail(plot_df$threshold, n=1)
  #print(t_end)
  
  step_inc <- 0.05   # x axis tick increment (change dynamically based on 
                     # range of threshold values eventually)
  
  if(v_count == 0){
    stop("Empty vertices in vertex.count. Ending plot utility")
  }
  
  factor <- e_count / v_count
  
  # ggplot handles colors by mapping the color to a string constant.
  # Thus, these variables won't work for keeping name changes consistent
  # edge_count_name <- "Edge Count"
  # vertex_count_name <- "Vertex Count"
  
  # Ensure that title is centered later
  ggplot2::theme_update(plot.title=ggplot2::element_text(hjust = 0.5))
  
 
  
  # ------------- OLD METHODS ------
  # unique_threshold <- plot_df$threshold[ !dup_vals ]
  # unique_edge_count <- plot_df$edge.count[ !dup_vals ]
  # unique_vertex_count <- plot_df$vertex.count [ !dup_vals ]
  # 
  # unique_df <- data.frame(threshold = unique_threshold,
  #                         edge.count = unique_edge_count,
  #                         vertex.count = unique_vertex_count)
  
  
  
  # Select all columns from df with unique threshold values (edge and vertex)
  # counts will also be the same because of the same iterative thresholding method
  dup_vals <- duplicated(plot_df[,c('threshold', 'edge.count', 'vertex.count')])
  unique_df <- plot_df [ !dup_vals , ]
  
  #print(unique_df)
  # ADD DOCUMENTATION TO EACH STEP!
  PLOT <- ggplot2::ggplot(data=unique_df,
                           ggplot2::aes(x=threshold)) + 
    
          ggplot2::geom_line(
            ggplot2::aes(y=vertex.count, 
                         color="Vertex Count")
          ) +
    
          ggplot2::geom_line(ggplot2::aes(y=edge.count / factor, 
                                          color="Edge Count")
                             ) + 
    

          ggplot2::xlab("Threshold Value") +
    
          ggplot2::ylab("Edge Count") +
    
          ggplot2::ggtitle("Edge and Vertex Count by Threshold Value") + 
    
          ggplot2::scale_x_continuous(
                          breaks = seq(t_begin, t_end, step_inc)
                        ) +
    
          ggplot2::scale_y_continuous(
                      name = "Vertex Count", 
                      sec.axis = ggplot2::sec_axis(~.*factor, name="Edge Count")
                    ) +
    
          ggplot2::scale_color_manual(name="Legend",
                                      breaks=c("Edge Count", "Vertex Count"),
                                      values=c("Edge Count" = "red",
                                               "Vertex Count" = "blue"
                                               )
                                      ) 
  
  # Add in markers for each method
  # print(names(labels))
  for(method_name in names(labels)){
    
    annot_string <- ""
    for(label_method in labels[[method_name]]){
      annot_string <- paste0(annot_string, label_method, '\n')
    }
    x_coord <- annotations[[method_name]]
    
    # Don't break the ggplot2::geom_point by a method having a suggested threshold
    # over the max possible threshold number. Not sure how this is happening for 
    # method='rmt' in particular, but need to double check get_results() is behaving
    # how it should.
    if(x_coord > t_end){
      next;
    }
    #print(paste0("edges: ", plot_df$edge.count[plot_df$threshold==x_coord]))
    # print(paste0("WE ARE NOW HERE: ", x_coord, method_name ))
    #print(paste0(method_name, "   ", unique_df$vertex.count[unique_df$threshold==x_coord]))
    #print(paste0(method_name, "   ", unique_df$edge.count[unique_df$threshold==x_coord]))
    
    PLOT <- PLOT +
            ggplot2::geom_vline(xintercept=x_coord,
                                linetype="solid",
                                color="orange",
                                linewidth=0.75
                                ) +

            # Point for the vertex curve
            # Shape 23 = 45 degree square
    
            ggplot2::geom_point(x=x_coord,
                                y=unique_df$vertex.count[unique_df$threshold==x_coord],
                                color="blue",
                                shape=23,
                                fill="blue",
                                size=3
                                ) +
            # # Point for the edge curve
            # Shape 23 = 45 degree square
            ggplot2::geom_point(x=x_coord,
                                y=(1/factor)*unique_df$edge.count[unique_df$threshold==x_coord],
                                color="red",
                                shape=23,
                                fill="red",
                                size=3
                                ) +
            # Add the method name annotations to the vertical line
            # Do this last so that the string is on the highest layer
            ggplot2::annotate("text",
                              label=annot_string,
                              x=x_coord,
                              y=0.1*v_count,
                              angle=25,
                              size=3.2
            )
  } # end of method annotation loop
  # print("done with loop")
  # Return the ggplot2 object for variable storage or direct plotting with show()
  return(PLOT)
  
  # Translated roughly from the following matplotlib code in Carissa's package (see DESCRIPTION file)
  #
  # number vertices and number edges vs thresholds
  # with sns.plotting_context("paper"):
  #   fig, ax = plt.subplots(figsize=(12, 6)) # long, high
  # 
  # ax.plot(df_plot["vertex-count"], alpha=0.6, label="Vertex count")
  # 
  # ax.plot(np.nan, np.nan, alpha=0.6, label="Edge-count", color="orange")
  # 
  # ax_twin = ax.twinx()
  # ax_twin.plot(df_plot["edge-count"], alpha=0.6, color="orange", label="Edge count")
  # 
  # xmin, xmax, ymin, ymax = ax.axis()
  # y_txt = (ymax - ymin)/100 + 10
  # 
  # for key, value in annotations.items():
  #   ax.plot(value, y_txt, 
  #           marker=11, #"CARETDOWNBASE"
  #           markersize=10,
  #           color="red")
  # 
  # txt = " &\n".join(labels[key])
  # ax.annotate(txt, (value, y_txt), 
  #             rotation=-45, 
  #             horizontalalignment='right', 
  #             verticalalignment='bottom')
  # 
  # ax.set(xlabel="Threshold")
  # ax.set(ylabel='Vertex count')
  # ax_twin.set(ylabel="Edge count")
  # ax_twin.grid(False)
  # 
  # ax.legend()
  # ax.legend(loc=1)
  # sns.despine(ax=ax,      right=True, left=False, bottom=True, top=True)
  # 
  # plt.tight_layout()
  
}  # end of plot_t_vs_ev()

