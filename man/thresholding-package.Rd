\name{thresholding}
\alias{thresholding-package}
\alias{thresholding}
\docType{package}
\title{
  Graph thresholding techniques and analysis methods for weighted graphs
}
\description{
  An R wrapper of C++ code that encapsulates several graph thresholding utilities and analysis tools. 
    Thresholding a graph helps to iteratively remove lower significance edges from a graph, and potentially shrinking the graph by pruning off disconnected vertices. 
    These methods readily apply to biological data, Pearson/Spearman correlations, and other types of graphs. 
    Our analysis tools take in a graph input file, and the user specifies what operations to perform on the graph to help determine the optimal threshold for the graph. 
    These thresholding tools rely on the igraph C library for graph creation and manipulation. Other functionalities are used from the alglib mathematical library as well. 
    After thresholding analysis is completed, a user can plot their results, extract the results from respective thresholding methods using get_results(), or the data frames produced from a specified method (iterative, local global, and significance values).
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Levi Hochstetler (lhochste@vols.utk.edu)

Maintainer: Levi Hochstetler (lhochste@vols.utk.edu)
}
\references{
  Bleker, Carissa & Grady, Stephen & Langston, Michael. (2023). A Brief Study of Gene Co-expression Thresholding Algorithms. 10.1007/978-981-99-7074-2_33. (https://www.researchgate.net/publication/374541767_A_Brief_Study_of_Gene_Co-expression_Thresholding_Algorithms)
  
  https://github.com/carissableker/thresholding
}
\keyword{ package }
\seealso{
  https://github.com/carissableker/thresholding
  https://github.com/snapplelevi/thresholding 
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
