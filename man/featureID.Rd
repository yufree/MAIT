\name{featureID}
\alias{featureID}
\title{Feature ID extractor from a MAIT object}
\description{
Function \code{featureID} extracts the feature IDs of a MAIT object as a vector. 
}
\usage{

featureID(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A numeric vector containing the feature IDs of the MAIT object.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
featureID(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
