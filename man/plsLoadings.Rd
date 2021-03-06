\name{plsLoadings}
\alias{plsLoadings}
\title{Loadings extractor for the PLS model}
\description{
Function \code{plsLoadings} returns the loading vectors for the PLS model when function \code{\link{plotPLS}} 
}
\usage{

plsLoadings(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}

}
\value{
A matrix with the PLS loading vectors.
}
\seealso{
  \code{\link{plotPLS}}}

  \examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPLS(MAIT,plot3d=FALSE)
plsLoadings(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}