\name{featureInfo}
\alias{featureInfo}
\title{Feature Info extractor from a MAIT object}
\description{
Function \code{featureInfo} extracts the slot MAIT.FeatureInfo of a MAIT object.
}
\usage{

featureInfo(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
An object of the class MAIT.FeatureInfo. More Info at \link{MAIT-class}.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
featureInfo(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
