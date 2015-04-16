\name{plsModel}
\alias{plsModel}
\title{Model extractor for either PLS}
\description{
Function \code{model} returns the model for PLS when function \code{\link{plotPLS}} have been already respectively launched. 
}
\usage{

plsModel(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
The PLS model of the MAIT.object
}
\seealso{
  \code{\link{plotPLS}}}

  \examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPLS(MAIT,plot3d=FALSE)
plsModel(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}