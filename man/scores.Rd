\name{scores}
\alias{scores}
\title{Retrieves the scores from a MAIT object}
\description{
Function \code{scores} extracts the scores MAIT object
}
\usage{

scores(object,
       type,...)

}

\arguments{
  \item{object}{
A \link{MAIT-class} object
}
  \item{type}{
If it is set to "none", the peak aggregation models are returned.
If it is set to "PCA", the PCA model is returned.
If it is set to "PLS", the PCA model is returned.

}
\item{...}{Other input}

}
\value{
A numeric matrix containing the scores saved in the MAIT object
}
\seealso{
  \code{\link{MAIT-class}} }

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
scores(MAIT)

MAIT<-plotPCA(MAIT,plot3d=FALSE)
scores(MAIT,type="PCA")

MAIT<-plotPLS(MAIT,plot3d=FALSE)
scores(MAIT,type="PLS")
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
