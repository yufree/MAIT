\name{pcaModel}
\alias{pcaModel}
\title{Model extractor for either PCA}
\description{
Function \code{model} returns the model for PCA when function \code{\link{plotPCA}} have been already respectively launched. 
}
\usage{

pcaModel(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
The PCA model of the MAIT.object
}
\seealso{
  \code{\link{plotPCA}}}
  
  \examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPCA(MAIT,plot3d=FALSE)
pcaModel(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}