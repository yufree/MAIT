\name{pcaScores}
\alias{pcaScores}
\title{Scores extractor for the PCA model}
\description{
Function \code{pcaScores} returns the loading vectors for the PCA model when function \code{\link{plotPCA}} 
}
\usage{

pcaScores(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}

}
\value{
A matrix with the PCA loading vectors.
}
\seealso{
  \code{\link{plotPCA}}}

  \examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPCA(MAIT,plot3d=FALSE)
pcaScores(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}