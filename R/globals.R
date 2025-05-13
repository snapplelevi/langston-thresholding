# Make sure that the column names from the output files from analysis()
# do not cause the 
#         "no visible binding for global variable XXXXXX"
# error in R CMD check...
#
# Recommended as a semi-hacky solution here: 
#     https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("threshold",            # Headers for analysis output files
                           "vertex.count",
                           "edge.count",
                           "connected.component.count",
                           "density",
                           "density.orig.V",
                           "largest.cc.size",
                           "second.largest.cc-size",
                           "clustering.coefficient",
                           "random.clustering.coefficient",
                           "second.eigenvalue",
                           "almost.disconnected.component.count",
                           "maximal.clique.count",
                           "clique.number",
                           "poisson.chi2",
                           "poisson.pvalue",
                           "goe.chi2",
                           "goe.pvalue",
                           "scale.free.KS",
                           "scale.free.KS.p.value",
                           "scale.free.alpha",     # End Headers for analysis output files
                           "maximal.clique.ratio", # Start/end dplyr variables in get_results
                           "threshold_column",     # R/absolute_threshold
                           "weight_column"         # R/edge_hist
                           ))
}


																	