% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummaryCKM.R
\name{summary.ckm}
\alias{summary.ckm}
\title{The function offers a summary of the results of CKM. The CKM analysis could either include the selection the number of clusters or without such an selection.}
\usage{
\method{summary}{ckm}(ckm.object)
}
\arguments{
\item{ckm.object}{an object of the S3 class "ckm", which could eitehr be the output of the function CKMSelVar (without the selection of number of clusters) or of the function CKMselAll (with the selection of nunber of clusters)}
}
\value{
a list of five items that are corresponding to the selected or determine number of clusters. The first and second are integers that indicate the number of clusters and masking variables, respectively.
The third is a vector that contains the indices of signaling variables. The forth demonstrates the cluster assignment. The last one, a matrix, represent the centroids of the signaling variables (columns) for all clusters (rows).
}
\description{
The function offers a summary of the results of CKM. The CKM analysis could either include the selection the number of clusters or without such an selection.
}
\examples{
ckm.sel.no <- CKMSelNo(dataset, ncluster, nnoisevar)
ckm.sel.var <- CKMSelVar(dataset, ncluster)
ckm.sel.all <- CKMSelAll(dataset, maxcluster)
summary(ckm.sel.no)
summary(ckm.sel.var)
summary(ckm.sel.all)


}
