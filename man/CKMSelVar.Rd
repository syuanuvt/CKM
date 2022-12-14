% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CKMSelVar.R
\name{CKMSelVar}
\alias{CKMSelVar}
\title{The function, combined with CardKMeans, selects the number of masking variables, given the fixed number of clusters}
\usage{
CKMSelVar(
  dataset,
  n.cluster,
  search = "dep",
  maxnum = 10,
  n.rep = 20,
  kmeans_starts = 10
)
}
\arguments{
\item{dataset}{the orginal dataset on which CKM and its model selection procedure operates}

\item{n.cluster}{the total number of clusters}

\item{search}{the mode of selecting over the grid. "all" = selecting over each point of the grid; while it maximizes the accuracy, it is overly slow with large number of variables.
"sub" = the "grid search with a zoom" strategy; while it is less accurate compared to searching the full grid, it is efficient even with large number of variables.
"dep" automatically adjust to one of the above two methods based on the number of variables. When # variables < 25, the search covers every possible value of the grid. This is also the default option.}

\item{maxnum}{the parameter is only useful when the "grid search with a zoom" strategy is applied. It restricts the maximal number of values searched over in any iteration. The default value is set at 10.}

\item{n.rep}{the number of permutated datasets when calculating the gap statistic}

\item{kmeans_starts}{the number of starts used in the kmeans algorithm}
}
\value{
@return The function will return a ckm object that is the list of five elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters; the fifth is the original dataset.
}
\description{
The function, combined with CardKMeans, selects the number of masking variables, given the fixed number of clusters
}
\examples{
ncluster <- 3
ckm.sel.var <- CKMSelVar(dataset, ncluster)

}
