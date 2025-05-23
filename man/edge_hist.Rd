% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edge_hist.R
\name{edge_hist}
\alias{edge_hist}
\title{edge_hist}
\usage{
edge_hist(infile, bin_width = 0.01, sep = "", title = "")
}
\arguments{
\item{infile}{string. File path to the .ncol formatted graph file (list of weighted edges - \strong{wel})}

\item{bin_width}{numeric. The width of each histogram bin. The default is bins of \strong{width} = 0.01}

\item{sep}{string. Specifically decide how the columns in the \strong{wel} file are separated.
The default separator is any white space.}

\item{title}{string. Optional user-supplied name for title. Default title is "Edge Edge Weight Histogram for <INPUT_FILE>"}
}
\value{
\code{ggplot2} object for modification and for showing with \code{show()}.
}
\description{
Histogram function for displaying edge weight frequencies. \code{edge_hist()} requires that an 
\url{https://lgl.sourceforge.net/} file be passed. In other words, the .ncol file can be thought
of as a \strong{w}eighted \strong{e}dge \strong{l}ist (\strong{wel}). The first column contains rows of vertices, the
second column contains rows of adjacent vertices, and the third column represents the weight for an edge
between the vertex in the first column one and the vertex in the second column.
}
\examples{
library(ThresholdTuner)
data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "ThresholdTuner") 
plot_obj <- edge_hist(data_file, bin_width = 0.01)
show(plot_obj)

}
