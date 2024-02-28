#' Plots the edge (E) and vertex (V)count of a graph at varying thresholds
#' Marks non-null methods and their optimal thresholds against the 
#' V/E line plots
#'
#' @param plot_df Dataframe with threshold values - returned from 
#' output of get_iterative_results, which contains the detailed
#' analysis of the graph at each increment of the thresholding
#' process.
#' @param D List of methods and their optimal thresholds from 
#' calling get_results. The user can either pass the resulting variable
#' from get_results or the list itself (i.e. D=variable$D instead of D=variable).
#' @export 
plot_t_vs_ev <- function(plot_df, D){
 
  annotations <- list()
  labels <- list()
  
  # Check if user passed the direct list from get_results() (D$D) or if they
  # passed just the outer variable D.
  if( ! ("ARTFUL.CHECK" %in% names(D)) ){
    D <- D$D
  }
  
  # Then remove ARTFUL.CHECK from the list before running the annotations
  # and labels loop
  D['ARTFUL.CHECK'] <- NULL
  
  index <- 1
  # Loop through items of D
  for(cur_val in D){
    method_name <- names(D)[index]
    print(paste0(method_name, ": ", cur_val))
    # If the method produced a non-zero value, add the method name (annotation)
    # and its corresponding value to the plotting sets.
    if(is.nan(cur_val) == FALSE && is.infinite(cur_val) == FALSE){

      is_in <- FALSE

      # Used to access annotation  name in annotations list
      annot_index <- 1

      for(annot_val in annotations){
        # String value of annotation name
        annot <- names(annotations)[annot_index]
        if(method_name == 'spectral_methods'){
          
          print(paste0(cur_val, " and ", annot_val))
          print(paste0("method_name: ", method_name, " and ", annot))
          print(dplyr::near(cur_val, annot_val, tol=1e-09))
        }        
        
        if(dplyr::near(cur_val, annot_val, tol=1e-09)){
          print(paste0("adding ", method_name, " to the ", annot, " vector"))
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
        print(paste0("adding: ", method_name))
        labels[[method_name]] <- c(method_name)   # start vector of labels, hashed on name
        annotations[[method_name]] <- cur_val    # value that corresponds to label name(s)
      }
      print(paste0("labels after", method_name))
      print(labels)
      print(paste0("annotations after ", method_name))
      print(annotations)
    }  # end outer is.nan if

    # Increment the index to get the string name of the method
    index <- index + 1

  } # end for loop
  
  
  print("Annotations: ")
  print(annotations)
  
  print("labels: ")
  print(labels)
  
  #print(head(plot_df))
  
  v_count <- plot_df$vertex.count[1]
  e_count <- plot_df$edge.count[1]
  t_begin <- plot_df$threshold[1]
  t_end <- utils::tail(plot_df$threshold, n=1)
  
  step_inc <- 0.05     # x axis tick increment
  
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
  
  # ADD DOCUMENTATION TO EACH STEP!
  PLOT <- ggplot2::ggplot(data=plot_df,
                           ggplot2::aes(x=threshold)) + 
          # Dummy comment
          ggplot2::geom_line(ggplot2::aes(y=edge.count / factor, 
                                          color="Edge Count")
                             ) + 
          ggplot2::geom_line(
            ggplot2::aes(y=vertex.count, 
                         color="Vertex Count")
            ) +
          ggplot2::xlab("Threshold Value") +
          ggplot2::ylab("Edge Count") +
          ggplot2::ggtitle("Edge and Vertex Count by Threshold Value") + 
          ggplot2::scale_x_continuous(
            breaks = seq(t_begin, t_end, step_inc)
          ) +
          ggplot2::scale_y_continuous(
            "Vertex Count", 
            sec.axis = ggplot2::sec_axis(~.*factor, name="Edge Count")
          ) +
          ggplot2::scale_color_manual(name="Legend",
                                      breaks=c("Edge Count", "Vertex Count"),
                                      values=c("Edge Count" = "red",
                                               "Vertex Count" = "blue"
                                               )
                                      ) 
  
  # Add in markers for each method
  for(method_name in names(labels)){
    
    annot_string <- ""
    for(label_method in labels[[method_name]]){
      annot_string <- paste0(annot_string, label_method, '\n')
    }
    x_coord <- annotations[[method_name]]
    #print(paste0("edges: ", plot_df$edge.count[plot_df$threshold==x_coord]))
    
    PLOT <- PLOT + 
            ggplot2::geom_vline(xintercept=x_coord,
                                linetype="solid",
                                color="orange",
                                linewidth=0.75
                                ) + 
      
            # Point for the vertex curve
            # Shape 23 = 45 degree square
            ggplot2::geom_point(x=x_coord,
                                y=plot_df$vertex.count[plot_df$threshold==x_coord],
                                color="blue",
                                shape=23,
                                fill="blue",
                                size=3
                                ) + 
            # Point for the edge curve
            # Shape 23 = 45 degree square
            ggplot2::geom_point(x=x_coord,
                                y=(1/factor)*plot_df$edge.count[plot_df$threshold==x_coord],
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
  
  # methods::show()
  show(PLOT)
  
  # Translated roughly from the following matplotlib code
  # Plotting stuff down here
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