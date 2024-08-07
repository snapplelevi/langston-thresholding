% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_results.R
\name{get_results}
\alias{get_results}
\title{Returns the resulting analysis method thresholding values
after running \code{thresholding::analysis()}}
\usage{
get_results(outfile_prefix, plot_iterative = FALSE)
}
\arguments{
\item{outfile_prefix}{File path and prefix for resulting output file from 
running the \code{analysis()} function (file would be <prefix>.iterative.txt). The 
output files will be assumed to be in the current working directory from where
\code{get_results()} was called if no path is specified.}

\item{plot_iterative}{Optionally plot the vertices and edges vs. threshold value as a graph.
uses ggplot2 to automatically call this package's plot_t_vs_ev() function without
the user needing to manually extract the required parameters.}
}
\value{
The returned list that can contain:
\itemize{
   \item A nested list of keyed on analysis method names. The values of these
   keys will be the recommended threshold from the corresponding method. This will 
   be called \code{"D"}.
   \item An alpha value representing the recommended significance value. This will 
   be called \code{"alpha"}
}
The output of \code{get_results()} and its separate wrappers will depend on the methods 
passed to \code{analysis()}. Values will either be valid or +/-\code{Inf}. If a method and/or value
for alpha does not show up in the returned list, then \code{get_results()} could not 
find this method or alpha value in the files with the matching prefix.
not show up
}
\description{
\strong{Note: modifying the column names of the output file(s) produced by} \code{analysis()} 
\strong{will lead to unintended behavior.}

\code{get_results()} will take an output file prefix and attempt to accumulate
the results connected to that prefix. This prefix come in the form of analysis methods
and their recommended thresholding value (as determined by \code{get_results()}).
}
\details{
For example, take the output file prefix of "\code{EXAMPLE_NAME}". If this is passed
to \code{analysis()} once with \code{methods=c(4,5)}, and another time with
\code{methods=c(6,7)}, then there will be two output files:
\enumerate{
   \item EXAMPLE_NAME-45.iterative.txt
   \item EXAMPLE_NAME-67.iterative.txt
}
When given the prefix of "\code{EXAMPLE_NAME}", \code{get_results()} will look at these 
files and combine the results from each method into one data structure.
}
\examples{
date_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
file.copy(data_file, "./")     # Copy the file to your working directory
outfile_prefix = "./get_results_test"
analysis(data_file, 
         outfile_prefix = outfile_prefix,
         methods = c(1,3,8),
         num_samples = 13,
         overwrite = TRUE
         )
thresholding::get_results(outfile_prefix, plot_iterative = TRUE)
}
