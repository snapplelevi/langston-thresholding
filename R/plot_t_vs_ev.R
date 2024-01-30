#' Plots the edge (E) and vertex (V)count of a graph at varying thresholds
#' Marks non-null methods and their optimal thresholds against the 
#' V/E line plots
#'
#' @param plot_df Dataframe with threshold values - returned from 
#' output of get_iterative_results, which contains the detailed
#' analysis of the graph at each increment of the thresholding
#' process
#' @param D List of methods and their optimal thresholds from 
#' calling get_results. The user can either pass the resulting variable
#' from get_results or the list itself (i.e. D=variable$D instead of D=variable).
#' @export 
plot_t_vs_ev <- function(plot_df, D){
  
  # UNCOMMENT OUT AFTER PLOT TESTING
  annotations <- list()
  labels <- list()

  index <- 1
  # Loop through items of D
  for(cur_val in D){
    method_name <- names(D)[index]

    # If the method produced a non-zero value, add the method name (annotation)
    # and its corresponding value to the plotting sets.
    if(is.nan(cur_val) == FALSE){

      is_in <- FALSE

      # Used to access annotation  name in annotations list
      annot_index <- 1

      for(annot_val in annotations){
        # String value of annotation name
        annot <- names(annotations)[annot_index]

        if(dplyr::near(cur_val, annot_val, tol=1e-09)){
          labels$annot <- append(annotations$annot, method_name)
          annotations$annot <- mean(cur_val, annot_val)

          # Make sure the element does not get added to the labels
          # since it is close to the current value and already in the set
          is_in <- TRUE

          # Increment name index as well
          annot_index <- annot_index + 1
        }
      }

      # Only add if the method hasn't been seen before
      if(is_in == FALSE){
        labels$method_name <- c(method_name)  # start list of labels, hashed on name
        annotations$method_name <- cur_val    # value that corresponds to label name(s)
      }

    }  # end outer is.nan if

    # Increment the index to get the string name of the method
    index <- index + 1

  } # end for loop
  print(head(plot_df))
  
  v_count <- plot_df$vertex.count[1]
  e_count <- plot_df$edge.count[1]
  
  if(v_count == 0){
    stop("Empty vertices in vertex.count. Ending plot utility")
  }
  
  factor <- e_count / v_count
  
  # ggplot handles colors by mapping the color to a string constant.
  # Thus, these variables won't work for keeping name changes consistent
  # edge_count_name <- "Edge Count"
  # vertex_count_name <- "Vertex Count"
  
  # ADD DOCUMENTATION TO EACH STEP!
  ggplot2::ggplot(data=plot_df,
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
  
}