\name{loadings}
\alias{loadings}
\title{Loadings extractor for either PCA or PLS models}
\description{
Function \code{loadings} returns the loading vectors for either the PCA, PLS models when functions \code{\link{plotPCA}} or \code{\link{plotPLS}} have been already respectively launched. It also can be used to retrieve the peak aggregation models.
}
\usage{

loadings(object, type = "none", ...)

}

\arguments{
  \item{object}{
A \link{MAIT-class} object
}
  \item{type}{
A character whose value should be "PCA" or "PLS" depending on which loading vectors are wanted. If it is set to "none", the peak aggregation models are retrieved.
}

\item{...}{Other input}

}
\value{
A matrix with the loading vectors.
}
\seealso{
  \code{\link{plotPCA}} or \code{\link{plotPLS}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
loadings(MAIT)

MAIT<-plotPCA(MAIT,plot3d=FALSE)
loadings(MAIT,type="PCA")

MAIT<-plotPLS(MAIT,plot3d=FALSE)
loadings(MAIT,type="PLS")
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
